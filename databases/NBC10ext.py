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

# <<< Supplementary Quantum Chemical Results >>>
DATA = {}

DATA['NUCLEAR REPULSION ENERGY'] = {}
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.2-dimer'          ] =     652.58240326
DATA['NUCLEAR REPULSION ENERGY']['NBC1-Bz-mono-unCP'              ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.3-dimer'          ] =     647.08083072
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.4-dimer'          ] =     641.79881504
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.5-dimer'          ] =     636.72435401
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.6-dimer'          ] =     631.84627841
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.7-dimer'          ] =     627.15417831
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.8-dimer'          ] =     622.63833806
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.9-dimer'          ] =     618.28967853
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.0-dimer'          ] =     614.09970566
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.1-dimer'          ] =     610.06046424
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.2-dimer'          ] =     606.16449631
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.5-dimer'          ] =     595.26834684
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-5.0-dimer'          ] =     579.39688238
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-5.5-dimer'          ] =     565.87021271
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-6.0-dimer'          ] =     554.22625379
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-6.5-dimer'          ] =     544.11253672
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-10.0-dimer'         ] =     499.16037479
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.4-dimer'          ] =     613.04854518
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.5-dimer'          ] =     608.81636557
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.6-dimer'          ] =     604.74550671
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.7-dimer'          ] =     600.82787505
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.8-dimer'          ] =     597.05577907
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.9-dimer'          ] =     593.42192782
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.0-dimer'          ] =     589.91942332
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.1-dimer'          ] =     586.54174882
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.2-dimer'          ] =     583.28275414
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.3-dimer'          ] =     580.13663931
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.4-dimer'          ] =     577.09793714
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.5-dimer'          ] =     574.16149552
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.6-dimer'          ] =     571.32245963
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.0-dimer'          ] =     560.85272572
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.5-dimer'          ] =     549.47925556
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.0-dimer'          ] =     539.65622514
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.5-dimer'          ] =     531.09189940
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-8.0-dimer'          ] =     523.56205991
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.2-dimer'       ] =     641.59153721
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.4-dimer'       ] =     640.97218086
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.6-dimer'       ] =     639.94808010
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.8-dimer'       ] =     638.53114770
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.0-dimer'       ] =     636.73745247
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.2-dimer'       ] =     634.58670201
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.4-dimer'       ] =     632.10168144
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.5-dimer'       ] =     630.74164257
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.6-dimer'       ] =     629.30768985
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.7-dimer'       ] =     627.80329032
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.8-dimer'       ] =     626.23200316
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.9-dimer'       ] =     624.59746513
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.0-dimer'       ] =     622.90337667
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.2-dimer'       ] =     619.35158842
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.4-dimer'       ] =     615.60701452
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.6-dimer'       ] =     611.70022314
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.8-dimer'       ] =     607.66157487
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-3.0-dimer'       ] =     603.52082284
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.2-dimer'           ] =     332.50866690
DATA['NUCLEAR REPULSION ENERGY']['NBC1-H2S-mono-unCP'             ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.4-dimer'           ] =     326.76493049
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.5-dimer'           ] =     324.08312886
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.6-dimer'           ] =     321.51823084
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.7-dimer'           ] =     319.06348175
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.8-dimer'           ] =     316.71257239
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.9-dimer'           ] =     314.45961051
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.0-dimer'           ] =     312.29909326
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.1-dimer'           ] =     310.22588084
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.2-dimer'           ] =     308.23517159
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.5-dimer'           ] =     302.71463310
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.75-dimer'          ] =     298.57449040
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.0-dimer'           ] =     294.79763877
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.25-dimer'          ] =     291.34045574
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.5-dimer'           ] =     288.16568982
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.0-dimer'           ] =     282.54011405
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.5-dimer'           ] =     277.71464354
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.0-dimer'           ] =     273.53417452
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.5-dimer'           ] =     269.88029141
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.2-dimer'            ] =     277.70122037
DATA['NUCLEAR REPULSION ENERGY']['NBC1-Bz2-mono-unCP'             ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-Me-mono-unCP'              ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.3-dimer'            ] =     276.14505886
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.4-dimer'            ] =     274.65657480
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.5-dimer'            ] =     273.23211647
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.6-dimer'            ] =     271.86820659
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.7-dimer'            ] =     270.56154682
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.8-dimer'            ] =     269.30901798
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.9-dimer'            ] =     268.10767718
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.0-dimer'            ] =     266.95475267
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.1-dimer'            ] =     265.84763738
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.2-dimer'            ] =     264.78388141
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.4-dimer'            ] =     262.77738579
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.6-dimer'            ] =     260.91850385
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.8-dimer'            ] =     259.19247204
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.0-dimer'            ] =     257.58628148
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.2-dimer'            ] =     256.08845607
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.4-dimer'            ] =     254.68885527
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.6-dimer'            ] =     253.37850109
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-6.0-dimer'            ] =     250.99455064
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.2-dimer'            ] =      42.94051671
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.3-dimer'            ] =      42.46449704
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.4-dimer'            ] =      42.01471911
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.5-dimer'            ] =      41.58914043
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.6-dimer'            ] =      41.18591734
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.7-dimer'            ] =      40.80338247
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.8-dimer'            ] =      40.44002498
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.9-dimer'            ] =      40.09447330
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.0-dimer'            ] =      39.76547998
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.1-dimer'            ] =      39.45190844
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.2-dimer'            ] =      39.15272123
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.3-dimer'            ] =      38.86696980
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.4-dimer'            ] =      38.59378540
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.6-dimer'            ] =      38.08199453
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.8-dimer'            ] =      37.61171219
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.0-dimer'            ] =      37.17815187
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.4-dimer'            ] =      36.40542136
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.8-dimer'            ] =      35.73746090
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.1-dimer'         ] =     664.74968142
DATA['NUCLEAR REPULSION ENERGY']['NBC1-Py-mono-unCP'              ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.3-dimer'         ] =     653.28897360
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.4-dimer'         ] =     647.90584891
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.5-dimer'         ] =     642.73711461
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.6-dimer'         ] =     637.77107423
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.7-dimer'         ] =     632.99683541
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.8-dimer'         ] =     628.40424073
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.9-dimer'         ] =     623.98380628
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.0-dimer'         ] =     619.72666684
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.1-dimer'         ] =     615.62452662
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.2-dimer'         ] =     611.66961499
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.3-dimer'         ] =     607.85464633
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.4-dimer'         ] =     604.17278378
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.5-dimer'         ] =     600.61760611
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.7-dimer'         ] =     593.86352067
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-5.0-dimer'         ] =     584.54275675
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-5.5-dimer'         ] =     570.86466240
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-6.0-dimer'         ] =     559.10620798
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-6.5-dimer'         ] =     548.90465922
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-7.0-dimer'         ] =     539.98032943
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.1-dimer'         ] =     631.74018099
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.3-dimer'         ] =     622.28221702
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.5-dimer'         ] =     613.57422251
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.6-dimer'         ] =     609.47520868
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.7-dimer'         ] =     605.53368830
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.8-dimer'         ] =     601.74111111
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.9-dimer'         ] =     598.08951503
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.0-dimer'         ] =     594.57147649
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.1-dimer'         ] =     591.18006603
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.2-dimer'         ] =     587.90880856
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.3-dimer'         ] =     584.75164753
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.4-dimer'         ] =     581.70291245
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.5-dimer'         ] =     578.75728949
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.7-dimer'         ] =     573.15574951
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.0-dimer'         ] =     565.41165299
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.5-dimer'         ] =     554.01089095
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-7.0-dimer'         ] =     544.16644693
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-8.0-dimer'         ] =     528.04095562
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-9.0-dimer'         ] =     515.40150653
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.2-dimer'       ] =     652.35026383
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.4-dimer'       ] =     651.65685475
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.6-dimer'       ] =     650.51106101
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.8-dimer'       ] =     648.92723975
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.0-dimer'       ] =     646.92462020
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.2-dimer'       ] =     644.52659143
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.4-dimer'       ] =     641.75995892
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.5-dimer'       ] =     640.24755050
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.6-dimer'       ] =     638.65423207
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.7-dimer'       ] =     636.98400901
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.8-dimer'       ] =     635.24097954
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.9-dimer'       ] =     633.42931896
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.0-dimer'       ] =     631.55326486
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.2-dimer'       ] =     627.62515488
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.4-dimer'       ] =     623.49127864
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.6-dimer'       ] =     619.18640729
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.8-dimer'       ] =     614.74502815
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-3.0-dimer'       ] =     610.20089775
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.2-dimer'       ] =     631.66053374
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.4-dimer'       ] =     631.10536715
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.6-dimer'       ] =     630.18691177
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.8-dimer'       ] =     628.91516711
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.0-dimer'       ] =     627.30369102
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.2-dimer'       ] =     625.36921338
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.4-dimer'       ] =     623.13120361
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.5-dimer'       ] =     621.90509666
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.6-dimer'       ] =     620.61142042
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.7-dimer'       ] =     619.25317914
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.8-dimer'       ] =     617.83346514
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.9-dimer'       ] =     616.35544587
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.0-dimer'       ] =     614.82235130
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.2-dimer'       ] =     611.60409513
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.4-dimer'       ] =     608.20532569
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.6-dimer'       ] =     604.65291019
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.8-dimer'       ] =     600.97358989
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-3.0-dimer'       ] =     597.19362514
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.2-dimer'          ] =     652.58240326
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.2-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.3-dimer'          ] =     647.08083072
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.3-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.4-dimer'          ] =     641.79881504
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.4-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.5-dimer'          ] =     636.72435401
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.6-dimer'          ] =     631.84627841
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.6-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.7-dimer'          ] =     627.15417831
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.7-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.8-dimer'          ] =     622.63833806
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.8-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.9-dimer'          ] =     618.28967853
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.9-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.0-dimer'          ] =     614.09970566
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.0-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.1-dimer'          ] =     610.06046424
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.1-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.2-dimer'          ] =     606.16449631
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.2-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.5-dimer'          ] =     595.26834684
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-5.0-dimer'          ] =     579.39688238
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-5.0-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-5.5-dimer'          ] =     565.87021271
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-5.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-6.0-dimer'          ] =     554.22625379
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-6.0-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-6.5-dimer'          ] =     544.11253672
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-6.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-10.0-dimer'         ] =     499.16037479
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-10.0-monoA-CP'      ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.4-dimer'          ] =     613.04854518
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.4-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.4-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.5-dimer'          ] =     608.81636557
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.5-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.6-dimer'          ] =     604.74550671
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.6-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.6-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.7-dimer'          ] =     600.82787505
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.7-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.7-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.8-dimer'          ] =     597.05577907
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.8-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.8-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.9-dimer'          ] =     593.42192782
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.9-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.9-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.0-dimer'          ] =     589.91942332
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.0-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.0-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.1-dimer'          ] =     586.54174882
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.1-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.1-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.2-dimer'          ] =     583.28275414
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.2-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.2-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.3-dimer'          ] =     580.13663931
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.3-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.3-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.4-dimer'          ] =     577.09793714
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.4-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.4-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.5-dimer'          ] =     574.16149552
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.5-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.6-dimer'          ] =     571.32245963
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.6-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.6-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.0-dimer'          ] =     560.85272572
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.0-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.0-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.5-dimer'          ] =     549.47925556
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.5-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.0-dimer'          ] =     539.65622514
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.0-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.0-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.5-dimer'          ] =     531.09189940
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.5-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-8.0-dimer'          ] =     523.56205991
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-8.0-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-8.0-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.2-dimer'       ] =     641.59153721
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.4-dimer'       ] =     640.97218086
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.6-dimer'       ] =     639.94808010
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.8-dimer'       ] =     638.53114770
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.0-dimer'       ] =     636.73745247
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.2-dimer'       ] =     634.58670201
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.4-dimer'       ] =     632.10168144
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.5-dimer'       ] =     630.74164257
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.5-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.6-dimer'       ] =     629.30768985
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.7-dimer'       ] =     627.80329032
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.7-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.8-dimer'       ] =     626.23200316
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.9-dimer'       ] =     624.59746513
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.9-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.0-dimer'       ] =     622.90337667
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.2-dimer'       ] =     619.35158842
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.4-dimer'       ] =     615.60701452
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.6-dimer'       ] =     611.70022314
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.8-dimer'       ] =     607.66157487
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-3.0-dimer'       ] =     603.52082284
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-3.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.2-dimer'           ] =     332.50866690
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.2-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.2-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.4-dimer'           ] =     326.76493049
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.4-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.4-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.5-dimer'           ] =     324.08312886
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.5-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.5-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.6-dimer'           ] =     321.51823084
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.6-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.6-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.7-dimer'           ] =     319.06348175
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.7-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.7-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.8-dimer'           ] =     316.71257239
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.8-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.8-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.9-dimer'           ] =     314.45961051
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.9-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.9-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.0-dimer'           ] =     312.29909326
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.0-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.0-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.1-dimer'           ] =     310.22588084
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.1-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.1-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.2-dimer'           ] =     308.23517159
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.2-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.2-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.5-dimer'           ] =     302.71463310
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.5-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.5-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.75-dimer'          ] =     298.57449040
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.75-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.75-monoB-CP'       ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.0-dimer'           ] =     294.79763877
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.0-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.0-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.25-dimer'          ] =     291.34045574
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.25-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.25-monoB-CP'       ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.5-dimer'           ] =     288.16568982
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.5-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.5-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.0-dimer'           ] =     282.54011405
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.0-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.0-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.5-dimer'           ] =     277.71464354
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.5-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.5-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.0-dimer'           ] =     273.53417452
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.0-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.0-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.5-dimer'           ] =     269.88029141
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.5-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.5-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.2-dimer'            ] =     277.70122037
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.2-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.2-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.3-dimer'            ] =     276.14505886
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.3-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.3-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.4-dimer'            ] =     274.65657480
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.4-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.4-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.5-dimer'            ] =     273.23211647
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.5-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.5-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.6-dimer'            ] =     271.86820659
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.6-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.6-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.7-dimer'            ] =     270.56154682
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.7-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.7-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.8-dimer'            ] =     269.30901798
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.8-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.8-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.9-dimer'            ] =     268.10767718
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.9-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.9-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.0-dimer'            ] =     266.95475267
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.0-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.0-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.1-dimer'            ] =     265.84763738
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.1-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.1-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.2-dimer'            ] =     264.78388141
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.2-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.2-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.4-dimer'            ] =     262.77738579
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.4-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.4-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.6-dimer'            ] =     260.91850385
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.6-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.6-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.8-dimer'            ] =     259.19247204
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.8-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.8-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.0-dimer'            ] =     257.58628148
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.0-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.0-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.2-dimer'            ] =     256.08845607
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.2-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.2-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.4-dimer'            ] =     254.68885527
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.4-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.4-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.6-dimer'            ] =     253.37850109
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.6-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.6-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-6.0-dimer'            ] =     250.99455064
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-6.0-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-6.0-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.2-dimer'            ] =      42.94051671
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.2-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.3-dimer'            ] =      42.46449704
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.3-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.4-dimer'            ] =      42.01471911
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.4-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.5-dimer'            ] =      41.58914043
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.5-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.6-dimer'            ] =      41.18591734
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.6-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.7-dimer'            ] =      40.80338247
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.7-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.8-dimer'            ] =      40.44002498
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.8-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.9-dimer'            ] =      40.09447330
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.9-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.0-dimer'            ] =      39.76547998
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.0-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.1-dimer'            ] =      39.45190844
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.1-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.2-dimer'            ] =      39.15272123
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.2-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.3-dimer'            ] =      38.86696980
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.3-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.4-dimer'            ] =      38.59378540
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.4-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.6-dimer'            ] =      38.08199453
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.6-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.8-dimer'            ] =      37.61171219
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.8-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.0-dimer'            ] =      37.17815187
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.0-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.4-dimer'            ] =      36.40542136
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.4-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.8-dimer'            ] =      35.73746090
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.8-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.1-dimer'         ] =     664.74968142
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.1-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.3-dimer'         ] =     653.28897360
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.3-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.4-dimer'         ] =     647.90584891
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.4-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.5-dimer'         ] =     642.73711461
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.5-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.6-dimer'         ] =     637.77107423
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.6-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.7-dimer'         ] =     632.99683541
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.7-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.8-dimer'         ] =     628.40424073
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.8-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.9-dimer'         ] =     623.98380628
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.9-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.0-dimer'         ] =     619.72666684
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.1-dimer'         ] =     615.62452662
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.1-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.2-dimer'         ] =     611.66961499
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.2-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.3-dimer'         ] =     607.85464633
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.3-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.4-dimer'         ] =     604.17278378
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.4-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.5-dimer'         ] =     600.61760611
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.5-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.7-dimer'         ] =     593.86352067
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.7-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-5.0-dimer'         ] =     584.54275675
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-5.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-5.5-dimer'         ] =     570.86466240
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-5.5-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-6.0-dimer'         ] =     559.10620798
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-6.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-6.5-dimer'         ] =     548.90465922
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-6.5-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-7.0-dimer'         ] =     539.98032943
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-7.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.1-dimer'         ] =     631.74018099
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.1-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.1-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.3-dimer'         ] =     622.28221702
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.3-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.3-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.5-dimer'         ] =     613.57422251
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.5-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.5-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.6-dimer'         ] =     609.47520868
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.6-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.6-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.7-dimer'         ] =     605.53368830
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.7-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.7-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.8-dimer'         ] =     601.74111111
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.8-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.8-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.9-dimer'         ] =     598.08951503
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.9-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.9-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.0-dimer'         ] =     594.57147649
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.0-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.1-dimer'         ] =     591.18006603
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.1-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.1-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.2-dimer'         ] =     587.90880856
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.2-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.2-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.3-dimer'         ] =     584.75164753
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.3-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.3-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.4-dimer'         ] =     581.70291245
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.4-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.4-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.5-dimer'         ] =     578.75728949
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.5-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.5-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.7-dimer'         ] =     573.15574951
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.7-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.7-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.0-dimer'         ] =     565.41165299
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.0-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.5-dimer'         ] =     554.01089095
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.5-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.5-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-7.0-dimer'         ] =     544.16644693
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-7.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-7.0-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-8.0-dimer'         ] =     528.04095562
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-8.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-8.0-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-9.0-dimer'         ] =     515.40150653
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-9.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-9.0-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.2-dimer'       ] =     652.35026383
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.4-dimer'       ] =     651.65685475
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.6-dimer'       ] =     650.51106101
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.8-dimer'       ] =     648.92723975
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.0-dimer'       ] =     646.92462020
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.2-dimer'       ] =     644.52659143
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.4-dimer'       ] =     641.75995892
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.5-dimer'       ] =     640.24755050
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.5-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.6-dimer'       ] =     638.65423207
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.7-dimer'       ] =     636.98400901
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.7-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.8-dimer'       ] =     635.24097954
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.9-dimer'       ] =     633.42931896
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.9-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.0-dimer'       ] =     631.55326486
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.2-dimer'       ] =     627.62515488
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.4-dimer'       ] =     623.49127864
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.6-dimer'       ] =     619.18640729
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.8-dimer'       ] =     614.74502815
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-3.0-dimer'       ] =     610.20089775
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-3.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.2-dimer'       ] =     631.66053374
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.4-dimer'       ] =     631.10536715
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.6-dimer'       ] =     630.18691177
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.8-dimer'       ] =     628.91516711
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.0-dimer'       ] =     627.30369102
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.2-dimer'       ] =     625.36921338
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.4-dimer'       ] =     623.13120361
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.5-dimer'       ] =     621.90509666
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.5-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.6-dimer'       ] =     620.61142042
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.7-dimer'       ] =     619.25317914
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.7-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.8-dimer'       ] =     617.83346514
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.9-dimer'       ] =     616.35544587
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.9-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.0-dimer'       ] =     614.82235130
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.2-dimer'       ] =     611.60409513
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.4-dimer'       ] =     608.20532569
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.6-dimer'       ] =     604.65291019
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.8-dimer'       ] =     600.97358989
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-3.0-dimer'       ] =     597.19362514
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-3.0-monoA-CP'    ] =     204.01997321
# DATA['SAPT MODELCHEM'] = 'SAPT3FC-SA-atz'
DATA['SAPT ELST ENERGY'] = {}
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-3.2'] =    -6.9579
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-3.3'] =    -4.7859
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-3.4'] =    -3.1971
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-3.5'] =    -2.0418
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-3.6'] =    -1.2082
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-3.7'] =    -0.6116
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-3.8'] =    -0.1893
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-3.9'] =     0.1052
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-4.0'] =     0.3066
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-4.1'] =     0.4401
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-4.2'] =     0.5246
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-4.5'] =     0.6042
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-5.0'] =     0.5178
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-5.5'] =     0.3941
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-6.0'] =     0.2935
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-6.5'] =     0.2191
DATA['SAPT ELST ENERGY']['NBC1-BzBz_S-10.0'] =     0.0381
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-4.3'] =    -9.0340
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-4.35'] =    -7.9365
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-4.4'] =    -6.9736
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-4.5'] =    -5.3910
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-4.6'] =    -4.1794
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-4.7'] =    -3.2544
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-4.8'] =    -2.5491
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-4.9'] =    -2.0117
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-5.0'] =    -1.6023
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-5.1'] =    -1.2900
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-5.2'] =    -1.0513
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-5.3'] =    -0.8682
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-5.4'] =    -0.7270
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-5.5'] =    -0.6175
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-5.6'] =    -0.5319
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-6.0'] =    -0.3298
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-6.5'] =    -0.2133
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-7.0'] =    -0.1489
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-7.5'] =    -0.1072
DATA['SAPT ELST ENERGY']['NBC1-BzBz_T-8.0'] =    -0.0786
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-0.2'] =    -3.1660
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-0.4'] =    -3.0799
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-0.6'] =    -2.9524
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-0.8'] =    -2.7994
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-1.0'] =    -2.6388
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-1.2'] =    -2.4863
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-1.4'] =    -2.3528
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-1.5'] =    -2.2944
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-1.6'] =    -2.2413
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-1.7'] =    -2.1930
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-1.8'] =    -2.1484
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-1.9'] =    -2.1064
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-2.0'] =    -2.0659
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-2.2'] =    -1.9834
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-2.4'] =    -1.8916
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-2.6'] =    -1.7844
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-2.8'] =    -1.6601
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD34-3.0'] =    -1.5213
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-3.15'] =   -10.3704
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-3.2'] =    -9.0920
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-3.3'] =    -7.0300
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-3.4'] =    -5.4653
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-3.5'] =    -4.2797
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-3.6'] =    -3.3812
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-3.7'] =    -2.6996
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-3.8'] =    -2.1812
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-3.9'] =    -1.7855
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-4.0'] =    -1.4817
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-4.1'] =    -1.2471
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-4.2'] =    -1.0645
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-4.5'] =    -0.7151
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-4.75'] =    -0.5513
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-5.0'] =    -0.4444
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-5.25'] =    -0.3687
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-5.5'] =    -0.3116
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-6.0'] =    -0.2299
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-6.5'] =    -0.1740
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-7.0'] =    -0.1340
DATA['SAPT ELST ENERGY']['NBC1-BzH2S-7.5'] =    -0.1045
DATA['SAPT ELST ENERGY']['NBC1-BzMe-3.15'] =    -4.7408
DATA['SAPT ELST ENERGY']['NBC1-BzMe-3.2'] =    -4.1123
DATA['SAPT ELST ENERGY']['NBC1-BzMe-3.3'] =    -3.1257
DATA['SAPT ELST ENERGY']['NBC1-BzMe-3.4'] =    -2.3753
DATA['SAPT ELST ENERGY']['NBC1-BzMe-3.5'] =    -1.8068
DATA['SAPT ELST ENERGY']['NBC1-BzMe-3.6'] =    -1.3771
DATA['SAPT ELST ENERGY']['NBC1-BzMe-3.7'] =    -1.0531
DATA['SAPT ELST ENERGY']['NBC1-BzMe-3.8'] =    -0.8093
DATA['SAPT ELST ENERGY']['NBC1-BzMe-3.9'] =    -0.6260
DATA['SAPT ELST ENERGY']['NBC1-BzMe-4.0'] =    -0.4882
DATA['SAPT ELST ENERGY']['NBC1-BzMe-4.1'] =    -0.3847
DATA['SAPT ELST ENERGY']['NBC1-BzMe-4.2'] =    -0.3068
DATA['SAPT ELST ENERGY']['NBC1-BzMe-4.4'] =    -0.2033
DATA['SAPT ELST ENERGY']['NBC1-BzMe-4.6'] =    -0.1429
DATA['SAPT ELST ENERGY']['NBC1-BzMe-4.8'] =    -0.1060
DATA['SAPT ELST ENERGY']['NBC1-BzMe-5.0'] =    -0.0822
DATA['SAPT ELST ENERGY']['NBC1-BzMe-5.2'] =    -0.0658
DATA['SAPT ELST ENERGY']['NBC1-BzMe-5.4'] =    -0.0540
DATA['SAPT ELST ENERGY']['NBC1-BzMe-5.6'] =    -0.0449
DATA['SAPT ELST ENERGY']['NBC1-BzMe-6.0'] =    -0.0320
DATA['SAPT ELST ENERGY']['NBC1-MeMe-3.1'] =    -1.4919
DATA['SAPT ELST ENERGY']['NBC1-MeMe-3.15'] =    -1.2506
DATA['SAPT ELST ENERGY']['NBC1-MeMe-3.2'] =    -1.0474
DATA['SAPT ELST ENERGY']['NBC1-MeMe-3.3'] =    -0.7324
DATA['SAPT ELST ENERGY']['NBC1-MeMe-3.4'] =    -0.5099
DATA['SAPT ELST ENERGY']['NBC1-MeMe-3.5'] =    -0.3532
DATA['SAPT ELST ENERGY']['NBC1-MeMe-3.6'] =    -0.2431
DATA['SAPT ELST ENERGY']['NBC1-MeMe-3.7'] =    -0.1661
DATA['SAPT ELST ENERGY']['NBC1-MeMe-3.8'] =    -0.1123
DATA['SAPT ELST ENERGY']['NBC1-MeMe-3.9'] =    -0.0749
DATA['SAPT ELST ENERGY']['NBC1-MeMe-4.0'] =    -0.0490
DATA['SAPT ELST ENERGY']['NBC1-MeMe-4.1'] =    -0.0312
DATA['SAPT ELST ENERGY']['NBC1-MeMe-4.2'] =    -0.0190
DATA['SAPT ELST ENERGY']['NBC1-MeMe-4.3'] =    -0.0108
DATA['SAPT ELST ENERGY']['NBC1-MeMe-4.4'] =    -0.0053
DATA['SAPT ELST ENERGY']['NBC1-MeMe-4.6'] =     0.0007
DATA['SAPT ELST ENERGY']['NBC1-MeMe-4.8'] =     0.0030
DATA['SAPT ELST ENERGY']['NBC1-MeMe-5.0'] =     0.0037
DATA['SAPT ELST ENERGY']['NBC1-MeMe-5.4'] =     0.0032
DATA['SAPT ELST ENERGY']['NBC1-MeMe-5.8'] =     0.0023
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-3.1'] =    -9.2779
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-3.2'] =    -6.7789
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-3.3'] =    -4.9427
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-3.4'] =    -3.5967
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-3.5'] =    -2.6126
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-3.6'] =    -1.8947
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-3.7'] =    -1.3723
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-3.8'] =    -0.9931
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-3.9'] =    -0.7188
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-4.0'] =    -0.5211
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-4.1'] =    -0.3791
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-4.2'] =    -0.2778
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-4.3'] =    -0.2059
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-4.4'] =    -0.1553
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-4.5'] =    -0.1200
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-4.7'] =    -0.0794
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-5.0'] =    -0.0578
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-5.5'] =    -0.0561
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-6.0'] =    -0.0610
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-6.5'] =    -0.0633
DATA['SAPT ELST ENERGY']['NBC1-PyPy_S2-7.0'] =    -0.0627
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-4.1'] =   -14.6558
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-4.27'] =    -9.4407
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-4.3'] =    -8.4738
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-4.36'] =    -7.2758
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-4.44'] =    -5.8622
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-4.5'] =    -5.0331
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-4.6'] =    -3.9380
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-4.7'] =    -3.1198
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-4.8'] =    -2.5065
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-4.9'] =    -2.0448
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-5.0'] =    -1.6951
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-5.1'] =    -1.4284
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-5.2'] =    -1.2231
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-5.3'] =    -1.0635
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-5.4'] =    -0.9379
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-5.5'] =    -0.8376
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-5.7'] =    -0.6894
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-6.0'] =    -0.5445
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-6.5'] =    -0.3986
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-7.0'] =    -0.3061
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-8.0'] =    -0.1933
DATA['SAPT ELST ENERGY']['NBC1-PyPy_T3-9.0'] =    -0.1299
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-0.2'] =    -6.8896
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-0.4'] =    -6.6980
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-0.6'] =    -6.4099
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-0.8'] =    -6.0602
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-1.0'] =    -5.6877
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-1.2'] =    -5.3271
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-1.4'] =    -5.0021
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-1.5'] =    -4.8559
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-1.6'] =    -4.7203
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-1.7'] =    -4.5941
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-1.8'] =    -4.4753
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-1.9'] =    -4.3613
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-2.0'] =    -4.2497
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-2.2'] =    -4.0229
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-2.4'] =    -3.7772
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-2.6'] =    -3.5021
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-2.8'] =    -3.1969
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD32-3.0'] =    -2.8684
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-0.2'] =    -1.1958
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-0.4'] =    -1.1626
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-0.6'] =    -1.1151
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-0.8'] =    -1.0607
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-1.0'] =    -1.0063
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-1.2'] =    -0.9589
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-1.4'] =    -0.9228
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-1.5'] =    -0.9095
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-1.6'] =    -0.8990
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-1.7'] =    -0.8912
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-1.8'] =    -0.8854
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-1.9'] =    -0.8810
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-2.0'] =    -0.8774
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-2.2'] =    -0.8698
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-2.4'] =    -0.8572
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-2.6'] =    -0.8362
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-2.8'] =    -0.8050
DATA['SAPT ELST ENERGY']['NBC1-BzBz_PD36-3.0'] =    -0.7644
DATA['SAPT EXCH ENERGY'] = {}
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-3.2'] =    24.9561
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-3.3'] =    18.7635
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-3.4'] =    14.0838
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-3.5'] =    10.5554
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-3.6'] =     7.8981
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-3.7'] =     5.9012
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-3.8'] =     4.4029
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-3.9'] =     3.2806
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-4.0'] =     2.4411
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-4.1'] =     1.8142
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-4.2'] =     1.3467
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-4.5'] =     0.5477
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-5.0'] =     0.1209
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-5.5'] =     0.0263
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-6.0'] =     0.0056
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-6.5'] =     0.0011
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_S-10.0'] =    -0.0000
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-4.3'] =    25.7550
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-4.35'] =    22.4627
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-4.4'] =    19.5678
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-4.5'] =    14.7981
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-4.6'] =    11.1425
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-4.7'] =     8.3556
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-4.8'] =     6.2415
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-4.9'] =     4.6452
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-5.0'] =     3.4453
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-5.1'] =     2.5469
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-5.2'] =     1.8769
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-5.3'] =     1.3792
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-5.4'] =     1.0107
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-5.5'] =     0.7387
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-5.6'] =     0.5387
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-6.0'] =     0.1490
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-6.5'] =     0.0286
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-7.0'] =     0.0052
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-7.5'] =     0.0009
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_T-8.0'] =     0.0002
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-0.2'] =    13.9153
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-0.4'] =    13.4339
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-0.6'] =    12.7058
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-0.8'] =    11.8244
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-1.0'] =    10.8888
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-1.2'] =     9.9837
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-1.4'] =     9.1649
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-1.5'] =     8.7968
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-1.6'] =     8.4568
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-1.7'] =     8.1432
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-1.8'] =     7.8529
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-1.9'] =     7.5818
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-2.0'] =     7.3252
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-2.2'] =     6.8352
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-2.4'] =     6.3463
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-2.6'] =     5.8313
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-2.8'] =     5.2777
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD34-3.0'] =     4.6880
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-3.15'] =    27.5031
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-3.2'] =    23.8335
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-3.3'] =    17.8287
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-3.4'] =    13.2948
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-3.5'] =     9.8849
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-3.6'] =     7.3294
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-3.7'] =     5.4209
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-3.8'] =     3.9997
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-3.9'] =     2.9446
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-4.0'] =     2.1633
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-4.1'] =     1.5862
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-4.2'] =     1.1609
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-4.5'] =     0.4506
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-4.75'] =     0.2027
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-5.0'] =     0.0905
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-5.25'] =     0.0401
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-5.5'] =     0.0177
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-6.0'] =     0.0034
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-6.5'] =     0.0006
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-7.0'] =     0.0001
DATA['SAPT EXCH ENERGY']['NBC1-BzH2S-7.5'] =     0.0000
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-3.15'] =    14.1126
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-3.2'] =    12.2598
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-3.3'] =     9.2408
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-3.4'] =     6.9365
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-3.5'] =     5.1863
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-3.6'] =     3.8634
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-3.7'] =     2.8678
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-3.8'] =     2.1216
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-3.9'] =     1.5647
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-4.0'] =     1.1505
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-4.1'] =     0.8436
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-4.2'] =     0.6170
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-4.4'] =     0.3276
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-4.6'] =     0.1724
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-4.8'] =     0.0899
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-5.0'] =     0.0464
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-5.2'] =     0.0238
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-5.4'] =     0.0121
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-5.6'] =     0.0061
DATA['SAPT EXCH ENERGY']['NBC1-BzMe-6.0'] =     0.0015
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-3.1'] =     5.1409
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-3.15'] =     4.3389
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-3.2'] =     3.6598
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-3.3'] =     2.5993
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-3.4'] =     1.8422
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-3.5'] =     1.3030
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-3.6'] =     0.9201
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-3.7'] =     0.6487
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-3.8'] =     0.4566
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-3.9'] =     0.3210
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-4.0'] =     0.2254
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-4.1'] =     0.1581
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-4.2'] =     0.1107
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-4.3'] =     0.0775
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-4.4'] =     0.0541
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-4.6'] =     0.0264
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-4.8'] =     0.0128
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-5.0'] =     0.0062
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-5.4'] =     0.0015
DATA['SAPT EXCH ENERGY']['NBC1-MeMe-5.8'] =     0.0003
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-3.1'] =    27.1989
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-3.2'] =    20.2791
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-3.3'] =    15.0925
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-3.4'] =    11.2132
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-3.5'] =     8.3174
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-3.6'] =     6.1600
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-3.7'] =     4.5554
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-3.8'] =     3.3640
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-3.9'] =     2.4808
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-4.0'] =     1.8271
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-4.1'] =     1.3440
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-4.2'] =     0.9875
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-4.3'] =     0.7248
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-4.4'] =     0.5316
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-4.5'] =     0.3895
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-4.7'] =     0.2088
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-5.0'] =     0.0817
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-5.5'] =     0.0169
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-6.0'] =     0.0034
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-6.5'] =     0.0007
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_S2-7.0'] =     0.0001
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-4.1'] =    42.6012
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-4.27'] =    25.4403
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-4.3'] =    23.2958
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-4.36'] =    19.4056
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-4.44'] =    15.1640
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-4.5'] =    12.5957
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-4.6'] =     9.2269
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-4.7'] =     6.7437
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-4.8'] =     4.9182
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-4.9'] =     3.5797
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-5.0'] =     2.6006
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-5.1'] =     1.8860
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-5.2'] =     1.3656
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-5.3'] =     0.9873
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-5.4'] =     0.7128
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-5.5'] =     0.5140
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-5.7'] =     0.2663
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-6.0'] =     0.0986
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-6.5'] =     0.0185
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-7.0'] =     0.0034
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-8.0'] =     0.0001
DATA['SAPT EXCH ENERGY']['NBC1-PyPy_T3-9.0'] =     0.0000
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-0.2'] =    24.6518
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-0.4'] =    23.7843
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-0.6'] =    22.4784
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-0.8'] =    20.9081
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-1.0'] =    19.2535
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-1.2'] =    17.6652
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-1.4'] =    16.2400
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-1.5'] =    15.6019
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-1.6'] =    15.0134
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-1.7'] =    14.4706
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-1.8'] =    13.9671
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-1.9'] =    13.4951
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-2.0'] =    13.0456
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-2.2'] =    12.1777
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-2.4'] =    11.2970
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-2.6'] =    10.3576
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-2.8'] =     9.3413
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD32-3.0'] =     8.2584
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-0.2'] =     7.8052
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-0.4'] =     7.5387
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-0.6'] =     7.1346
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-0.8'] =     6.6426
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-1.0'] =     6.1169
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-1.2'] =     5.6046
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-1.4'] =     5.1383
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-1.5'] =     4.9278
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-1.6'] =     4.7331
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-1.7'] =     4.5536
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-1.8'] =     4.3876
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-1.9'] =     4.2332
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-2.0'] =     4.0878
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-2.2'] =     3.8135
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-2.4'] =     3.5442
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-2.6'] =     3.2641
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-2.8'] =     2.9647
DATA['SAPT EXCH ENERGY']['NBC1-BzBz_PD36-3.0'] =     2.6460
DATA['SAPT IND ENERGY'] = {}
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-3.2'] =    -0.7208
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-3.3'] =    -0.6112
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-3.4'] =    -0.5156
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-3.5'] =    -0.4342
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-3.6'] =    -0.3659
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-3.7'] =    -0.3090
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-3.8'] =    -0.2619
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-3.9'] =    -0.2228
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-4.0'] =    -0.1904
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-4.1'] =    -0.1633
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-4.2'] =    -0.1406
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-4.5'] =    -0.0919
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-5.0'] =    -0.0481
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-5.5'] =    -0.0268
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-6.0'] =    -0.0155
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-6.5'] =    -0.0090
DATA['SAPT IND ENERGY']['NBC1-BzBz_S-10.0'] =    -0.0002
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-4.3'] =    -2.9732
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-4.35'] =    -2.6429
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-4.4'] =    -2.3470
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-4.5'] =    -1.8470
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-4.6'] =    -1.4511
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-4.7'] =    -1.1396
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-4.8'] =    -0.8957
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-4.9'] =    -0.7054
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-5.0'] =    -0.5569
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-5.1'] =    -0.4412
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-5.2'] =    -0.3509
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-5.3'] =    -0.2803
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-5.4'] =    -0.2249
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-5.5'] =    -0.1814
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-5.6'] =    -0.1469
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-6.0'] =    -0.0664
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-6.5'] =    -0.0285
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-7.0'] =    -0.0143
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-7.5'] =    -0.0080
DATA['SAPT IND ENERGY']['NBC1-BzBz_T-8.0'] =    -0.0047
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-0.2'] =    -0.5412
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-0.4'] =    -0.6077
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-0.6'] =    -0.6960
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-0.8'] =    -0.7849
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-1.0'] =    -0.8560
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-1.2'] =    -0.8981
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-1.4'] =    -0.9087
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-1.5'] =    -0.9033
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-1.6'] =    -0.8919
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-1.7'] =    -0.8756
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-1.8'] =    -0.8555
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-1.9'] =    -0.8324
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-2.0'] =    -0.8073
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-2.2'] =    -0.7533
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-2.4'] =    -0.6969
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-2.6'] =    -0.6391
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-2.8'] =    -0.5796
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD34-3.0'] =    -0.5180
DATA['SAPT IND ENERGY']['NBC1-BzH2S-3.15'] =    -4.2717
DATA['SAPT IND ENERGY']['NBC1-BzH2S-3.2'] =    -3.7447
DATA['SAPT IND ENERGY']['NBC1-BzH2S-3.3'] =    -2.8759
DATA['SAPT IND ENERGY']['NBC1-BzH2S-3.4'] =    -2.2084
DATA['SAPT IND ENERGY']['NBC1-BzH2S-3.5'] =    -1.6969
DATA['SAPT IND ENERGY']['NBC1-BzH2S-3.6'] =    -1.3056
DATA['SAPT IND ENERGY']['NBC1-BzH2S-3.7'] =    -1.0065
DATA['SAPT IND ENERGY']['NBC1-BzH2S-3.8'] =    -0.7782
DATA['SAPT IND ENERGY']['NBC1-BzH2S-3.9'] =    -0.6039
DATA['SAPT IND ENERGY']['NBC1-BzH2S-4.0'] =    -0.4708
DATA['SAPT IND ENERGY']['NBC1-BzH2S-4.1'] =    -0.3692
DATA['SAPT IND ENERGY']['NBC1-BzH2S-4.2'] =    -0.2914
DATA['SAPT IND ENERGY']['NBC1-BzH2S-4.5'] =    -0.1502
DATA['SAPT IND ENERGY']['NBC1-BzH2S-4.75'] =    -0.0916
DATA['SAPT IND ENERGY']['NBC1-BzH2S-5.0'] =    -0.0586
DATA['SAPT IND ENERGY']['NBC1-BzH2S-5.25'] =    -0.0391
DATA['SAPT IND ENERGY']['NBC1-BzH2S-5.5'] =    -0.0269
DATA['SAPT IND ENERGY']['NBC1-BzH2S-6.0'] =    -0.0137
DATA['SAPT IND ENERGY']['NBC1-BzH2S-6.5'] =    -0.0074
DATA['SAPT IND ENERGY']['NBC1-BzH2S-7.0'] =    -0.0043
DATA['SAPT IND ENERGY']['NBC1-BzH2S-7.5'] =    -0.0026
DATA['SAPT IND ENERGY']['NBC1-BzMe-3.15'] =    -1.3593
DATA['SAPT IND ENERGY']['NBC1-BzMe-3.2'] =    -1.2026
DATA['SAPT IND ENERGY']['NBC1-BzMe-3.3'] =    -0.9393
DATA['SAPT IND ENERGY']['NBC1-BzMe-3.4'] =    -0.7326
DATA['SAPT IND ENERGY']['NBC1-BzMe-3.5'] =    -0.5712
DATA['SAPT IND ENERGY']['NBC1-BzMe-3.6'] =    -0.4459
DATA['SAPT IND ENERGY']['NBC1-BzMe-3.7'] =    -0.3488
DATA['SAPT IND ENERGY']['NBC1-BzMe-3.8'] =    -0.2738
DATA['SAPT IND ENERGY']['NBC1-BzMe-3.9'] =    -0.2157
DATA['SAPT IND ENERGY']['NBC1-BzMe-4.0'] =    -0.1708
DATA['SAPT IND ENERGY']['NBC1-BzMe-4.1'] =    -0.1359
DATA['SAPT IND ENERGY']['NBC1-BzMe-4.2'] =    -0.1087
DATA['SAPT IND ENERGY']['NBC1-BzMe-4.4'] =    -0.0708
DATA['SAPT IND ENERGY']['NBC1-BzMe-4.6'] =    -0.0472
DATA['SAPT IND ENERGY']['NBC1-BzMe-4.8'] =    -0.0324
DATA['SAPT IND ENERGY']['NBC1-BzMe-5.0'] =    -0.0229
DATA['SAPT IND ENERGY']['NBC1-BzMe-5.2'] =    -0.0167
DATA['SAPT IND ENERGY']['NBC1-BzMe-5.4'] =    -0.0125
DATA['SAPT IND ENERGY']['NBC1-BzMe-5.6'] =    -0.0096
DATA['SAPT IND ENERGY']['NBC1-BzMe-6.0'] =    -0.0060
DATA['SAPT IND ENERGY']['NBC1-MeMe-3.1'] =    -0.2675
DATA['SAPT IND ENERGY']['NBC1-MeMe-3.15'] =    -0.2251
DATA['SAPT IND ENERGY']['NBC1-MeMe-3.2'] =    -0.1892
DATA['SAPT IND ENERGY']['NBC1-MeMe-3.3'] =    -0.1332
DATA['SAPT IND ENERGY']['NBC1-MeMe-3.4'] =    -0.0934
DATA['SAPT IND ENERGY']['NBC1-MeMe-3.5'] =    -0.0653
DATA['SAPT IND ENERGY']['NBC1-MeMe-3.6'] =    -0.0456
DATA['SAPT IND ENERGY']['NBC1-MeMe-3.7'] =    -0.0319
DATA['SAPT IND ENERGY']['NBC1-MeMe-3.8'] =    -0.0223
DATA['SAPT IND ENERGY']['NBC1-MeMe-3.9'] =    -0.0156
DATA['SAPT IND ENERGY']['NBC1-MeMe-4.0'] =    -0.0109
DATA['SAPT IND ENERGY']['NBC1-MeMe-4.1'] =    -0.0077
DATA['SAPT IND ENERGY']['NBC1-MeMe-4.2'] =    -0.0054
DATA['SAPT IND ENERGY']['NBC1-MeMe-4.3'] =    -0.0038
DATA['SAPT IND ENERGY']['NBC1-MeMe-4.4'] =    -0.0027
DATA['SAPT IND ENERGY']['NBC1-MeMe-4.6'] =    -0.0014
DATA['SAPT IND ENERGY']['NBC1-MeMe-4.8'] =    -0.0007
DATA['SAPT IND ENERGY']['NBC1-MeMe-5.0'] =    -0.0004
DATA['SAPT IND ENERGY']['NBC1-MeMe-5.4'] =    -0.0002
DATA['SAPT IND ENERGY']['NBC1-MeMe-5.8'] =    -0.0001
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-3.1'] =    -1.1221
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-3.2'] =    -0.8956
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-3.3'] =    -0.7160
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-3.4'] =    -0.5744
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-3.5'] =    -0.4630
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-3.6'] =    -0.3755
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-3.7'] =    -0.3065
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-3.8'] =    -0.2521
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-3.9'] =    -0.2089
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-4.0'] =    -0.1745
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-4.1'] =    -0.1469
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-4.2'] =    -0.1247
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-4.3'] =    -0.1066
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-4.4'] =    -0.0917
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-4.5'] =    -0.0795
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-4.7'] =    -0.0609
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-5.0'] =    -0.0424
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-5.5'] =    -0.0250
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-6.0'] =    -0.0154
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-6.5'] =    -0.0097
DATA['SAPT IND ENERGY']['NBC1-PyPy_S2-7.0'] =    -0.0063
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-4.1'] =    -4.8762
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-4.27'] =    -2.9633
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-4.3'] =    -2.7165
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-4.36'] =    -2.2840
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-4.44'] =    -1.8167
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-4.5'] =    -1.5322
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-4.6'] =    -1.1574
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-4.7'] =    -0.8786
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-4.8'] =    -0.6707
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-4.9'] =    -0.5154
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-5.0'] =    -0.3990
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-5.1'] =    -0.3114
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-5.2'] =    -0.2452
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-5.3'] =    -0.1948
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-5.4'] =    -0.1562
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-5.5'] =    -0.1265
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-5.7'] =    -0.0851
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-6.0'] =    -0.0502
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-6.5'] =    -0.0243
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-7.0'] =    -0.0135
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-8.0'] =    -0.0052
DATA['SAPT IND ENERGY']['NBC1-PyPy_T3-9.0'] =    -0.0023
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-0.2'] =    -0.7740
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-0.4'] =    -0.9133
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-0.6'] =    -1.0970
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-0.8'] =    -1.2784
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-1.0'] =    -1.4192
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-1.2'] =    -1.4981
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-1.4'] =    -1.5133
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-1.5'] =    -1.5003
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-1.6'] =    -1.4764
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-1.7'] =    -1.4439
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-1.8'] =    -1.4050
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-1.9'] =    -1.3619
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-2.0'] =    -1.3162
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-2.2'] =    -1.2215
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-2.4'] =    -1.1260
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-2.6'] =    -1.0304
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-2.8'] =    -0.9326
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD32-3.0'] =    -0.8313
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-0.2'] =    -0.3783
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-0.4'] =    -0.4105
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-0.6'] =    -0.4533
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-0.8'] =    -0.4970
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-1.0'] =    -0.5327
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-1.2'] =    -0.5547
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-1.4'] =    -0.5611
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-1.5'] =    -0.5587
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-1.6'] =    -0.5530
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-1.7'] =    -0.5446
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-1.8'] =    -0.5338
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-1.9'] =    -0.5212
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-2.0'] =    -0.5072
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-2.2'] =    -0.4760
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-2.4'] =    -0.4424
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-2.6'] =    -0.4070
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-2.8'] =    -0.3702
DATA['SAPT IND ENERGY']['NBC1-BzBz_PD36-3.0'] =    -0.3321
DATA['SAPT DISP ENERGY'] = {}
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-3.2'] =   -13.2507
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-3.3'] =   -11.3615
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-3.4'] =    -9.7658
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-3.5'] =    -8.4149
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-3.6'] =    -7.2670
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-3.7'] =    -6.2897
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-3.8'] =    -5.4554
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-3.9'] =    -4.7417
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-4.0'] =    -4.1295
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-4.1'] =    -3.6035
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-4.2'] =    -3.1503
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-4.5'] =    -2.1302
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-5.0'] =    -1.1541
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-5.5'] =    -0.6566
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-6.0'] =    -0.3913
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-6.5'] =    -0.2431
DATA['SAPT DISP ENERGY']['NBC1-BzBz_S-10.0'] =    -0.0190
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-4.3'] =   -11.8637
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-4.35'] =   -10.9763
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-4.4'] =   -10.1588
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-4.5'] =    -8.7106
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-4.6'] =    -7.4789
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-4.7'] =    -6.4295
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-4.8'] =    -5.5346
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-4.9'] =    -4.7705
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-5.0'] =    -4.1179
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-5.1'] =    -3.5599
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-5.2'] =    -3.0825
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-5.3'] =    -2.6739
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-5.4'] =    -2.3239
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-5.5'] =    -2.0236
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-5.6'] =    -1.7660
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-6.0'] =    -1.0471
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-6.5'] =    -0.5739
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-7.0'] =    -0.3329
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-7.5'] =    -0.2032
DATA['SAPT DISP ENERGY']['NBC1-BzBz_T-8.0'] =    -0.1294
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-0.2'] =    -9.7308
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-0.4'] =    -9.6230
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-0.6'] =    -9.4419
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-0.8'] =    -9.1949
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-1.0'] =    -8.8917
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-1.2'] =    -8.5430
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-1.4'] =    -8.1589
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-1.5'] =    -7.9568
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-1.6'] =    -7.7493
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-1.7'] =    -7.5375
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-1.8'] =    -7.3220
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-1.9'] =    -7.1036
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-2.0'] =    -6.8828
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-2.2'] =    -6.4361
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-2.4'] =    -5.9848
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-2.6'] =    -5.5324
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-2.8'] =    -5.0818
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD34-3.0'] =    -4.6374
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-3.15'] =   -30.5311
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-3.2'] =   -10.4662
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-3.3'] =    -8.8654
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-3.4'] =    -7.5246
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-3.5'] =    -6.3984
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-3.6'] =    -5.4505
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-3.7'] =    -4.6514
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-3.8'] =    -3.9769
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-3.9'] =    -3.4063
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-4.0'] =    -2.9232
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-4.1'] =    -2.5133
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-4.2'] =    -2.1652
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-4.5'] =    -1.4018
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-4.75'] =    -0.9907
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-5.0'] =    -0.7101
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-5.25'] =    -0.5165
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-5.5'] =    -0.3812
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-6.0'] =    -0.2165
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-6.5'] =    -0.1294
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-7.0'] =    -0.0807
DATA['SAPT DISP ENERGY']['NBC1-BzH2S-7.5'] =    -0.0522
DATA['SAPT DISP ENERGY']['NBC1-BzMe-3.15'] =  -119.9394
DATA['SAPT DISP ENERGY']['NBC1-BzMe-3.2'] =    -6.4603
DATA['SAPT DISP ENERGY']['NBC1-BzMe-3.3'] =    -5.5160
DATA['SAPT DISP ENERGY']['NBC1-BzMe-3.4'] =    -4.7150
DATA['SAPT DISP ENERGY']['NBC1-BzMe-3.5'] =    -4.0347
DATA['SAPT DISP ENERGY']['NBC1-BzMe-3.6'] =    -3.4565
DATA['SAPT DISP ENERGY']['NBC1-BzMe-3.7'] =    -2.9647
DATA['SAPT DISP ENERGY']['NBC1-BzMe-3.8'] =    -2.5463
DATA['SAPT DISP ENERGY']['NBC1-BzMe-3.9'] =    -2.1899
DATA['SAPT DISP ENERGY']['NBC1-BzMe-4.0'] =    -1.8863
DATA['SAPT DISP ENERGY']['NBC1-BzMe-4.1'] =    -1.6276
DATA['SAPT DISP ENERGY']['NBC1-BzMe-4.2'] =    -1.4069
DATA['SAPT DISP ENERGY']['NBC1-BzMe-4.4'] =    -1.0570
DATA['SAPT DISP ENERGY']['NBC1-BzMe-4.6'] =    -0.8005
DATA['SAPT DISP ENERGY']['NBC1-BzMe-4.8'] =    -0.6115
DATA['SAPT DISP ENERGY']['NBC1-BzMe-5.0'] =    -0.4712
DATA['SAPT DISP ENERGY']['NBC1-BzMe-5.2'] =    -0.3664
DATA['SAPT DISP ENERGY']['NBC1-BzMe-5.4'] =    -0.2875
DATA['SAPT DISP ENERGY']['NBC1-BzMe-5.6'] =    -0.2276
DATA['SAPT DISP ENERGY']['NBC1-BzMe-6.0'] =    -0.1463
DATA['SAPT DISP ENERGY']['NBC1-MeMe-3.1'] =    -2.9665
DATA['SAPT DISP ENERGY']['NBC1-MeMe-3.15'] =    -2.7005
DATA['SAPT DISP ENERGY']['NBC1-MeMe-3.2'] =    -2.4599
DATA['SAPT DISP ENERGY']['NBC1-MeMe-3.3'] =    -2.0448
DATA['SAPT DISP ENERGY']['NBC1-MeMe-3.4'] =    -1.7039
DATA['SAPT DISP ENERGY']['NBC1-MeMe-3.5'] =    -1.4235
DATA['SAPT DISP ENERGY']['NBC1-MeMe-3.6'] =    -1.1924
DATA['SAPT DISP ENERGY']['NBC1-MeMe-3.7'] =    -1.0017
DATA['SAPT DISP ENERGY']['NBC1-MeMe-3.8'] =    -0.8439
DATA['SAPT DISP ENERGY']['NBC1-MeMe-3.9'] =    -0.7131
DATA['SAPT DISP ENERGY']['NBC1-MeMe-4.0'] =    -0.6045
DATA['SAPT DISP ENERGY']['NBC1-MeMe-4.1'] =    -0.5141
DATA['SAPT DISP ENERGY']['NBC1-MeMe-4.2'] =    -0.4386
DATA['SAPT DISP ENERGY']['NBC1-MeMe-4.3'] =    -0.3755
DATA['SAPT DISP ENERGY']['NBC1-MeMe-4.4'] =    -0.3226
DATA['SAPT DISP ENERGY']['NBC1-MeMe-4.6'] =    -0.2404
DATA['SAPT DISP ENERGY']['NBC1-MeMe-4.8'] =    -0.1815
DATA['SAPT DISP ENERGY']['NBC1-MeMe-5.0'] =    -0.1387
DATA['SAPT DISP ENERGY']['NBC1-MeMe-5.4'] =    -0.0838
DATA['SAPT DISP ENERGY']['NBC1-MeMe-5.8'] =    -0.0528
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-3.1'] =   -14.1088
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-3.2'] =   -12.0373
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-3.3'] =   -10.2965
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-3.4'] =    -8.8280
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-3.5'] =    -7.5866
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-3.6'] =    -6.5349
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-3.7'] =    -5.6412
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-3.8'] =    -4.8806
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-3.9'] =    -4.2315
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-4.0'] =    -3.6768
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-4.1'] =    -3.2016
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-4.2'] =    -2.7935
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-4.3'] =    -2.4426
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-4.4'] =    -2.1405
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-4.5'] =    -1.8798
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-4.7'] =    -1.4591
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-5.0'] =    -1.0139
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-5.5'] =    -0.5760
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-6.0'] =    -0.3433
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-6.5'] =    -0.2135
DATA['SAPT DISP ENERGY']['NBC1-PyPy_S2-7.0'] =    -0.1377
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-4.1'] =   -15.1126
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-4.27'] =  -558.9631
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-4.3'] =   -10.7116
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-4.36'] =  -234.2433
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-4.44'] =    -8.4847
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-4.5'] =    -7.6911
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-4.6'] =    -6.5439
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-4.7'] =    -5.5816
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-4.8'] =    -4.7720
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-4.9'] =    -4.0894
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-5.0'] =    -3.5127
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-5.1'] =    -3.0246
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-5.2'] =    -2.6107
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-5.3'] =    -2.2591
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-5.4'] =    -1.9597
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-5.5'] =    -1.7044
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-5.7'] =    -1.2994
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-6.0'] =    -0.8822
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-6.5'] =    -0.4870
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-7.0'] =    -0.2850
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-8.0'] =    -0.1125
DATA['SAPT DISP ENERGY']['NBC1-PyPy_T3-9.0'] =    -0.0511
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-0.2'] =   -13.1990
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-0.4'] =   -13.0383
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-0.6'] =   -12.7706
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-0.8'] =   -12.4107
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-1.0'] =   -11.9749
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-1.2'] =   -11.4806
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-1.4'] =   -10.9446
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-1.5'] =   -10.6650
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-1.6'] =   -10.3794
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-1.7'] =   -10.0890
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-1.8'] =    -9.7946
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-1.9'] =    -9.4968
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-2.0'] =    -9.1963
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-2.2'] =    -8.5885
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-2.4'] =    -7.9745
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-2.6'] =    -7.3576
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-2.8'] =    -6.7417
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD32-3.0'] =    -6.1320
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-0.2'] =    -7.2423
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-0.4'] =    -7.1666
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-0.6'] =    -7.0402
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-0.8'] =    -6.8670
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-1.0'] =    -6.6521
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-1.2'] =    -6.4026
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-1.4'] =    -6.1255
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-1.5'] =    -5.9788
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-1.6'] =    -5.8277
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-1.7'] =    -5.6728
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-1.8'] =    -5.5150
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-1.9'] =    -5.3546
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-2.0'] =    -5.1923
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-2.2'] =    -4.8633
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-2.4'] =    -4.5305
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-2.6'] =    -4.1974
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-2.8'] =    -3.8661
DATA['SAPT DISP ENERGY']['NBC1-BzBz_PD36-3.0'] =    -3.5396
