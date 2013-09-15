#import os
#import re
#import math
#import copy
from molecule import Molecule
#from periodictable import *
#from physconst import *
from vecutil import *
from exceptions import *
#from libmintscoordentry import *

#LINEAR_A_TOL = 1.0E-2  # When sin(a) is below this, we consider the angle to be linear
#DEFAULT_SYM_TOL = 1.0E-8
#FULL_PG_TOL = 1.0e-8
#ZERO = 1.0E-14
NOISY_ZERO = 1.0E-8
COORD_ZERO = 1.0E-5  # tolerance in coordinate alignment btwn qc programs


class OrientMols(object):
    """Class to encode a transformation between two molecular coordinate
    systems. After initializing with two qcdb.Molecule objects at the
    same geometry in possible different frames and orderings, class
    can apply the appropriate transformations to coordinate, gradient,
    Hessian, etc. arrays.

    """

    def __init__(self, molPermanent, molChangeable):
        """Stores the shift, rotation, axis exchange, axis inversion,
        and atom remapping necessary to bring the geometry of
        *molChangeable* into coincidence with the geometry of
        *molPermanent*. *molPermanent* and *molChangeable* must be
        :py:class:`qcdb.Molecule` and represent the same geometry.

        """
        # <<< Permanent (Psi4) >>>

        # Molecule
        self.Pmol = molPermanent
        # Vector to shift Pmol to center of mass
        self.Pshift = []
        # Matrix to rotate Pmol to inertial frame
        self.Protate = []

        # <<< Changeable (Cfour) >>>

        # Molecule
        self.Cmol = molChangeable
        # Vector to shift Cmol to center of mass
        self.Cshift = []
        # Matrix to rotate Cmol to inertial frame
        self.Crotate = []
        # Matrix to rotate Cmol to axis representation of Pmol
        self.Cexchflip = []
        # Vector to map Cmol to atom ordering of Pmol
        self.Catommap = []

        try:
            if ((self.Pmol.natom() == self.Cmol.natom()) and \
               (abs(self.Pmol.nuclear_repulsion_energy() - self.Cmol.nuclear_repulsion_energy()) < 1.0e-3) and \
               (self.Pmol.rotor_type() == self.Cmol.rotor_type())):
                self.create_orientation_from_molecules(self.Pmol, self.Cmol)
            else:
                print 'qcdb.orient.__init__ debug info'
                self.Pmol.print_out()
                print('natom', self.Pmol.natom(), 'NRE', self.Pmol.nuclear_repulsion_energy(), 'rotor', self.Pmol.rotor_type())
                self.Cmol.print_out()
                print('natom', self.Cmol.natom(), 'NRE', self.Cmol.nuclear_repulsion_energy(), 'rotor', self.Cmol.rotor_type())
                raise ValidationError("""OrientMols Molecule arguments differ fatally.""")
        except AttributeError:
            raise ValidationError("""OrientMols must be instantiated with two qcdb.Molecule objects.""")

    def __str__(self):
        text = """  ==> qcdb OrientMols <==\n\n"""
        text += """   natom:     %d\n\n""" % (self.Pmol.natom())
        text += """   PNRE:       %16.8f\n""" % (self.Pmol.nuclear_repulsion_energy())
        text += """   Pshift:    %s\n""" % (self.Pshift)
        text += """   Protate:   %s\n""" % (self.Protate)
        text += """\n   CNRE:       %16.8f\n""" % (self.Cmol.nuclear_repulsion_energy())
        text += """   Cshift:    %s\n""" % (self.Cshift)
        text += """   Crotate:   %s\n""" % (self.Crotate)
        text += """   Cexchflip: %s\n""" % (self.Cexchflip)
        text += """   Catommap:  %s\n""" % (self.Catommap)
        return text

    def create_orientation_from_molecules(self, Pmol, Cmol):
        """Finds the shift, rotation, axis exchange, axis inversion,
        and atom remapping necessary to bring the geometry of *Cmol*
        into coincidence with the geometry of *Pmol*. *Pmol* and *Cmol*
        must be :py:class:`qcdb.Molecule` and represent the same
        geometry. Presently catches some errors of orientation that
        Cfour as *Cmol* should properly fulfill. These are unnecessary
        restrictions and can be relaxed later.

        """
        p4mol = Pmol.clone()
        c4mol = Cmol.clone()
        Nat = p4mol.natom()
        eye3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

        # Find translation to CoM, straightforward
        com = p4mol.center_of_mass()
#        com = [ x if abs(x) > NOISY_ZERO else 0.0 for x in com]
        p4mol.translate(scale(com, -1.0))
        self.Pshift = com

        com = c4mol.center_of_mass()
#        com = [ x if abs(x) > NOISY_ZERO else 0.0 for x in com]
        c4mol.translate(scale(com, -1.0))
        self.Cshift = com
        #   Extra check since Cfour always returns at center of mass
        if not (all([abs(com[i]) < COORD_ZERO for i in range(3)])):
            print 'qcdb.orient.create_orientation_from_molecules debug info'
            print '\ncom', com
            raise ValidationError("""molChangeable not at center of mass.""")

        # Find rotation to MoI frame, straightforward
        moi, frame = p4mol.inertial_system(zero=NOISY_ZERO)
        Psort = sorted(range(3), key=lambda x: moi[x])
#        frame = [[ x if abs(x) > NOISY_ZERO else 0.0 for x in ax] for ax in frame]
        p4mol.rotate(frame)
        self.Protate = frame

        moi, frame = c4mol.inertial_system(zero=NOISY_ZERO)
        #moi, frame = c4mol.inertial_system(zero=1.0e-6)
        Csort = sorted(range(3), key=lambda x: moi[x])
#        frame = [[ x if abs(x) > NOISY_ZERO else 0.0 for x in ax] for ax in frame]
        c4mol.rotate(frame)
        self.Crotate = frame
        #   Extra check since Cfour always returns in inertial frame
        if not all([all([abs(frame[i][j] - eye3[i][j]) < COORD_ZERO for j in range(3)]) for i in range(3)]):
            print 'qcdb.orient.create_orientation_from_molecules debug info'
            print '\nCgeom: '
            for item in c4mol.geometry():
                print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
            print '\nmoi', moi, '\nframe', frame, '\nCsort', Csort
            raise ValidationError("""molChangeable not in inertial frame.""")

        # Find degrees of freedom among axis exchanges
        rotor = p4mol.rotor_type()
        if rotor == 'RT_ATOM':
            freebytop = []
        elif rotor == 'RT_LINEAR':                 # 0  <  IB == IC      inf > B == C
            freebytop = [1, 2]
        elif rotor == 'RT_SPHERICAL_TOP':          # IA == IB == IC       A == B == C
            freebytop = [0, 1, 2]
        elif rotor == 'RT_PROLATE_SYMMETRIC_TOP':  # IA <  IB == IC       A >  B == C
            freebytop = [1, 2]
        elif rotor == 'RT_OBLATE_SYMMETRIC_TOP':   # IA == IB <  IC       A == B >  C
            freebytop = [0, 1]
        elif rotor == 'RT_ASYMMETRIC_TOP':         # IA <  IB <  IC       A >  B >  C
            freebytop = []

        # Find mapping of axis exchange and flipping that brings Cgeom into coincidence with Pgeom
        Pgeom = p4mol.geometry()
        Cgeom = c4mol.geometry()

        exchMat = zero(3, 3)
        Pwhite = list(range(3))
        Cwhite = list(range(3))
        freeaxphase = [] #[1, 1, 1]]

        while len(Pwhite) > 0:
            Paxs = Pwhite[0]
            allowed = list(set(freebytop) & set(Cwhite)) if Paxs in freebytop else [Paxs]

            for Caxs in allowed:
                exchEntry = None
                PcolS = sorted([row[Psort[Paxs]] for row in Pgeom])
                CcolS = sorted([row[Csort[Caxs]] for row in Cgeom])
                CcolMS = sorted([-row[Csort[Caxs]] for row in Cgeom])

#                print [row[Psort[Paxs]] for row in Pgeom]
#                print [row[Csort[Caxs]] for row in Cgeom]
#                print PcolS
#                print CcolS
                if all([abs(CcolS[at] - PcolS[at]) < COORD_ZERO for at in range(Nat)]):
                    if all([abs(CcolS[at] - CcolMS[at]) < COORD_ZERO for at in range(Nat)]):
#                        PcolI = sorted(range(Nat), key=lambda x: [row[Psort[Paxs]] for row in Pgeom][x])
#                        CcolI = sorted(range(Nat), key=lambda x: [row[Csort[Caxs]] for row in Cgeom][x])
#                        print PcolI[0] > PcolI[-1], PcolI[0], PcolI[-1]
#                        print CcolI[0] > CcolI[-1], CcolI[0], CcolI[-1]
#                        if (PcolI[0] > PcolI[-1]) == (CcolI[0] > CcolI[-1]):
#                            exchEntry = 1  # Pos when symm P4 col == C4 col and pos vals at same col ends
#                            print 'a', exchEntry
#                        else:
#                            exchEntry = -1  # Neg when symm P4 col == C4 col and pos vals at diff col ends
#                            print 'b', exchEntry
                        exchEntry = 1  # Indeterminate when symm P4 col == C4 col
#                        cfap = freeaxphase[:]
#                        print freeaxphase, cfap
#                        for indx in rang
#                        for ph in cfap:
#                            ph[Psort[Paxs]] = -1
#                            freeaxphase.append(ph)
                        freeaxphase.append(Psort[Paxs])
                        print 'b'
                        print PcolS
                    else:
                        exchEntry = 1  # Pos when asym P4 col == C4 col
                        print 'c', exchEntry
                        print PcolS
                elif all([abs(CcolMS[at] - PcolS[at]) < COORD_ZERO for at in range(Nat)]):
                    exchEntry = -1  # Neg when P4 col == -C4 col
                    print 'd', exchEntry
                    print PcolS

                if exchEntry is not None:
                    exchMat[Csort[Caxs]][Psort[Paxs]] = exchEntry
                    Pwhite.remove(Paxs)
                    Cwhite.remove(Caxs)
                    break
            else:
                print 'qcdb.orient.create_orientation_from_molecules debug info'
                print '\nrotor', rotor, 'Paxs', Paxs, 'Caxs', Caxs, 'Pwhite', Pwhite, 'Cwhite', Cwhite, \
                    'allowed', allowed, 'P(axs)', Psort.index(Paxs), 'C(Caxs)', Csort.index(Caxs), \
                    '\nPcolS: ', PcolS, '\nCcolS: ', CcolS, '\nCcolMS: ', CcolMS
                print '\nexchMat'
                show(exchMat)
                print '\nPgeom: '
                for item in Pgeom:
                    print '       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2])
                print '\nCgeom: '
                for item in Cgeom:
                    print '       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2])
                print self
                raise ValidationError("""Axis unreconcilable between QC programs.""")

        self.Cexchflip = exchMat
        c4mol.rotate(exchMat)
        print 'freeaxphase', freeaxphase

#        ph2 = [[1, 1, 1]]
#        for ax in freeaxphase:
#            for indx in range(len(ph2)):
#                asdf = copy.deepcopy(ph2[indx])
#                asdf[ax] = -1
#                ph2.append(asdf)
#
#        print 'ph2', ph2
        show(exchMat)

        # Find mapping of atom exchange that brings Cgeom into coincidence with Pgeom
        Cgeom = c4mol.geometry()

#        print '\nPgeom: '
#        for item in Pgeom:
#            print '       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2])
#        print '\nCgeom: '
#        for item in Cgeom:
#            print '       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2])

        mapMat = [None] * Nat
        Pwhite = list(range(Nat))
        Cwhite = list(range(Nat))

        while len(Pwhite) > 0:
            Patm = Pwhite[0]
            sameElem = [at for at in range(Nat) if c4mol.symbol(at) == p4mol.symbol(Patm)]
            allowed = list(set(sameElem) & set(Cwhite))


            allowedPh = [[1, 1, 1]]
            for ax in freeaxphase:
                for indx in range(len(allowedPh)):
                    asdf = copy.deepcopy(allowedPh[indx])
                    asdf[ax] = -1
                    allowedPh.append(asdf)
                allowedPh = sorted(allowedPh, reverse=True, key=lambda x: sum(x))
#            print 'Patm', Patm, 'allowedPH', allowedPh

#            print '\nAAAAA Patm', Patm, 'Pwhite', Pwhite, \
#                'Cwhite', Cwhite, 'sameElem', sameElem, 'allowed', allowed, \
#                'allowedPh', allowedPh, \
#                '\nPgeom[Patm]: ', Pgeom[Patm], '\nmapMat', mapMat
            for ph in allowedPh:
                for Catm in allowed:
#                    print 'ph Catm', ph, Catm
                    print 'Pgeom', Pgeom[Patm]
                    print 'Cgeom', Cgeom[Catm]
                    print 'CgeomP', [Cgeom[Catm][ax] * ph[ax] for ax in range(3)]

                    if all([abs(Cgeom[Catm][ax] * ph[ax] - Pgeom[Patm][ax]) < COORD_ZERO for ax in range(3)]):

                        mapMat[Patm] = Catm
                        print 'setting mapMat[', Patm, '] = ', Catm, 'under phase', ph
                        for ax in range(3):
                            #if ph[ax] == -1 and abs(Pgeom[Patm][ax]) > COORD_ZERO:
                            if ax in freeaxphase and abs(Pgeom[Patm][ax]) > COORD_ZERO:
                                freeaxphase.remove(ax)
                                print 'removing', ax, 'now freeaxphase', freeaxphase
#                                print 'freeaxphaseREV', freeaxphase
                                for at in range(Nat):
                                    Cgeom[at][ax] *= ph[ax]
                                for ax2 in range(3):
                                    exchMat[ax2][ax] *= ph[ax]
                                print 'new Cexchflip', exchMat
                        Pwhite.remove(Patm)
                        Cwhite.remove(Catm)
                        break
                if mapMat[Patm] is not None:
                    break
            else:
                print 'qcdb.orient.create_orientation_from_molecules debug info'
                print '\nPatm', Patm, 'Catm', Catm, 'ph', ph, 'Pwhite', Pwhite, \
                    'Cwhite', Cwhite, 'sameElem', sameElem, 'allowed', allowed, \
                    'allowedPh', allowedPh, '\nCgeom[Catm]: ', Cgeom[Catm], \
                    '\nPgeom[Patm]: ', Pgeom[Patm], '\nmapMat', mapMat
                print '\nPgeom: '
                for item in Pgeom:
                    print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
                print('\nCgeom: ')
                for item in Cgeom:
                    print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
                print self
                raise ValidationError("""Atom unreconcilable between QC programs.""")

        self.Catommap = mapMat
        self.Cexchflip = exchMat
        # Note that this is resetting the geom but not the atoms, so c4mol.print_out() is deceptive
        new_geom = []
        for at in range(Nat):
            new_geom.append(Cgeom[mapMat[at]])
        c4mol.set_geometry(new_geom)

        # One last check that p4mol and c4mol align
        Pgeom = p4mol.geometry()
        Cgeom = c4mol.geometry()

        if not all([all([abs(Cgeom[at][ax] - Pgeom[at][ax]) < COORD_ZERO for ax in range(3)]) for at in range(Nat)]):
            raise ValidationError("""Geometries unreconcilable between QC programs:\n  P4 %s\n  C4 %s""" % (Pgeom, Cgeom))

    def transform_coordinates(self, coord):
        """

        """
        #print self
        print "Original"
        coord.print_out()

        coord.translate(scale(self.Cshift, -1.0))
        print "Shift"
        coord.print_out()

        coord.rotate(self.Crotate)
        print "Rotate"
        coord.print_out()

        coord.rotate(self.Cexchflip)
        print "ExchFlip"
        coord.print_out()

        geom = coord.geometry()
        new_geom = []
        for at in range(coord.natom()):
            new_geom.append(geom[self.Catommap[at]])
        coord.set_geometry(new_geom)
        print "AtomMap"
        coord.print_out()

        coord.rotate(transpose(self.Protate))
        print "P4 Rotate"
        coord.print_out()

        coord.translate(self.Pshift)
        print "P4 Shift"
        coord.print_out()

    def transform_coordinates2(self, coord):
        """

        """
        geom = coord.geometry()

#        print self
#        print "Original"
#        for item in geom:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        # ?????
        #coord.translate(scale(self.Cshift, -1.0))
        geom2 = []
        for item in geom:
            geom2.append(sub(item, self.Cshift))
        geom = geom2
#        print "Shift"
#        for item in geom:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        #coord.rotate(self.Crotate)
        geom2 = mult(geom, self.Crotate)
        geom = geom2
#        print "Rotate"
#        for item in geom:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        #coord.rotate(self.Cexchflip)
        geom2 = mult(geom, self.Cexchflip)
        geom = geom2
 #       print "ExchFlip"
 #       for item in geom:
 #           print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        geom2 = []
        for at in range(coord.natom()):
            geom2.append(geom[self.Catommap[at]])
        coord.set_geometry(geom2)
        geom = geom2
#        print "AtomMap"
#        for item in geom:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        #coord.rotate(transpose(self.Protate))
        geom2 = mult(geom, transpose(self.Protate))
        geom = geom2
#        print "P4 Rotate"
#        for item in geom:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        geom2 = []
        for item in geom:
            geom2.append(add(item, self.Pshift))
        geom = geom2
#        print "P4 Shift"
#        for item in geom:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        Pgeom = self.Pmol.geometry()
        Cgeom = geom
        Nat = len(geom)
        if not all([all([abs(Cgeom[at][ax] - Pgeom[at][ax]) < COORD_ZERO for ax in range(3)]) for at in range(Nat)]):
            raise ValidationError("""Geometries unreconcilable between QC programs:\n  P4 %s\n  C4 %s""" % (Pgeom, Cgeom))
        return geom

    def transform_gradient(self, arr):
        """Applies to *arr* the transformation appropriate to bring a
        gradient in *molChangeable* orientation into *molPermanent*
        orientation. In particular, applies a rotation to place it
        in the inertial frame, a column exchange and phasing to place
        it in the axis system, a row exchange to place it in the atom
        ordering, and a rotation to remove it from the inertial frame.

        """
        arr = mult(arr, self.Crotate)
#        print "Rotate"
#        for item in arr:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        arr = mult(arr, self.Cexchflip)
#        print "ExchFlip"
#        for item in arr:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        arr2 = []
        for at in range(len(arr)):
            arr2.append(arr[self.Catommap[at]])
        arr = arr2
#        print "AtomMap"
#        for item in arr:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        arr = mult(arr, transpose(self.Protate))
#        print "P4 Rotate"
#        for item in arr:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        return arr
