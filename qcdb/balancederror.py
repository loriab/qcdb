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

from __future__ import print_function

def compute_balanced_error(mcsys, refsys, refeq, curveratio, verbose=False):
    """Computes the balanced error quantity from the interaction energy
    of a system under investigation with a certain model chemistry,
    *mcsys*. Also required are certain figures about the dissociation
    curve computed with a reference model chemistry: the interaction
    energy of the system *refsys*, the interaction energy of the
    equilibrium system *refeq*, and a measure of the position of the
    system along the dissociation curve, generally the ratio of system
    to equilibrium separations *curveratio* (e.g., curveratio = 0.846
    for BzBz_S-3.3 (equilibrium is 3.9 A)). Notation follows suppmat
    of JCP ???. If *refeq* equals *refsys* and *curveratio* equals 1,
    returned quantity equals simple relative error.

    """
    if mcsys is None or refsys is None:
        return [None, None, None]

    # the following two parameters define the balanced relative error
    m = 0.03  # minimum permitted weight for a point
    p = 10.0  # multiples of |REFeq| above REFeq to which zero-line in head is displaced

    # compute q and weight from reference quantities
    if curveratio < 1.0:
        q = p  # head of curve
    else:
        q = 1.0  # equilibrium or tail of curve

    temp = q - 1.0 + refsys / refeq
    weight = max(m, temp / q)

    # compute error and weighted error from modelchem quantities
    error = (mcsys - refsys) * q / abs(refeq * temp)
    xsys = weight * error
    if verbose:
        print("""%65s wtdamp %8.3f   posndamp %8.3f   damp %8.3f %8.3f""" % ('',
            xsys / ((mcsys - refsys) / abs(refeq)),
            abs(refsys) / abs(refeq),
            (xsys / ((mcsys - refsys) / abs(refeq))) * (abs(refsys)/abs(refeq)),
            xsys / ((mcsys - refsys)/abs(refsys))))

    # [absolute, relative, balanced error]
    return [mcsys - refsys, (mcsys - refsys) / abs(refsys), xsys]


if __name__ == '__main__':

    # BzBz_S B3LYP-D3/aTZ
    refeq = -1.72
    refR = 3.9
    data = [
        [3.1, 40.0, 42.0],  # this data point is fake
        [3.15, 20.0, 22.0],  # ditto
        [3.2, 3.462, 4.41],
        [3.3, 1.484, 2.05],
        [3.4, 0.147, 0.48],
        [3.5, -0.724, -0.53],
        [3.6, -1.259, -1.14],
        [3.7, -1.558, -1.48],
        [3.8, -1.693, -1.65],
        [3.9, -1.717, -1.69],
        [4.0, -1.693, -1.66],
        [4.1, -1.577, -1.58],
        [4.2, -1.459, -1.48],
        [4.5, -1.066, -1.10],
        [5.0, -0.546, -0.58],
        [5.5, -0.251, -0.30],
        [6.0, -0.101, -0.16],
        [6.5, -0.029, -0.09],
        [10.0, 0.018, 0.01]]

    print('%10s %10s %10s %10s %10s %10s' %
        ('R', 'refIE', 'mcIE', 'Error', '%Err', '%BalErr'))
    for sys in data:
        ber = compute_balanced_error(sys[2], sys[1], refeq, sys[0] / refR, verbose=True)
        print('%10.1f %10.2f %10.2f %10.2f %10.1f %10.1f' % \
            (sys[0], sys[1], sys[2], ber[0], 100 * ber[1], 100 * ber[2]))
