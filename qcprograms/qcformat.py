"""Parent classes for quantum chemistry program input and output file
formats.
"""
import re


class InputFormat(object):

    def __init__(self, mem, mtd, bas, mol, sys, cast):

        # total job memory in MB
        self.memory = mem
        # computational method
        self.method = mtd.lower()
        # qcdb.Molecule object
        self.molecule = mol
        # database member index
        self.index = sys
        # orbital basis set
        self.basis = bas.lower()
        # do cast up from sto-3g basis?
        self.castup = cast

    def corresponding_aux_basis(self):
        """For Dunning basis sets, returns strings from which auxiliary
        basis sets and heavy-aug can be constructed. Note that
        valence/core-valence/etc. is conserved and X-zeta/(X+d)zeta is
        not, since this is the usual aux basis pattern.
        *augbasis* is round up to the nearest aug-cc-pVXZ
        *rootbasis* is round down to the nearest cc-pVXZ
        *auxbasis* is round up to the nearest cc-pVXZ or aug-cc-pVXZ
        """
        Dunmatch = re.compile(r'^(.*cc-)(pv|pcv|pwcv).*?([dtq56]).*z$').match(self.basis)

        if Dunmatch:
            rootbas = 'cc-' + Dunmatch.group(2) + Dunmatch.group(3) + 'z'
            augbas = 'aug-cc-' + Dunmatch.group(2) + Dunmatch.group(3) + 'z'
            if Dunmatch.group(1) == 'cc-':
                auxbas = rootbas
            else:
                auxbas = augbas
        else:
            rootbas = None
            augbas = None
            auxbas = None

        return [rootbas, augbas, auxbas]
