"""

"""
import re
import periodictable
from xcpt import *

atom = re.compile(r"""^(?P<isotope>[0-9]{1,3})?  # optional isotope
                       (?P<symbol>[A-Z]{1,3}) +  # mandatory atomic symbol
                       (?P<label>(_\w+)|(\d+))?$ # optional userid""",
                  re.IGNORECASE | re.X)


class AtomGist(object):
    """

    """

    def __init__(self, label, x, y, z, mass=None, charge=None, tags={}):
        """

        """
        #: label for position
        self.label = label

        #: element label for position (read-only, dep. on label)

        #: Atomic Number of position (read-only, dep. on label)

        #: x position in Bohr
        self.x = float(x)

        #: y position in Bohr
        self.y = float(y)

        #: z position in Bohr
        self.z = float(z)

        #: Nuclear charge on location in electron units
        self.charge = charge

        #: mass of position in Daltons
        self.mass = mass

        #: all tags associated with position
        self.tags = {}
        for item, val in tags.iteritems():
            self.tags[item.lower()] = val

    def __str__(self):
        return """  %-23s %16.8f %16.8f %16.8f    %6.2f    %8.3f %s""" % (
            self.label + ' / ' + self.symbol + ' / ' + str(self.Z),
            self.x, self.y, self.z, self.charge, self.mass,
            ', '.join(['{}={}'.format(k, v) for k, v in self.tags.iteritems()]))

    @property
    def label(self):
        """label for position"""
        ans = self.__iso if self.__iso else ''
        ans += self.symbol
        ans += self.__lbl if self.__lbl else ''
        return ans

    @label.setter
    def label(self, val):
        atm = atom.match(val.strip().capitalize())
        self.__iso = atm.group('isotope')
        self.__el = atm.group('symbol').upper()
        self.__lbl = atm.group('label')
        print '==>', self.__iso, self.__el, self.__lbl, '<=='
        if not atm:
            raise ValidationError(
                """Position label not of Co, Co3, 60Co (isotope), or """
                """Co_III form, case insensitive: {}.""".format(
                    val))
        if self.__el not in periodictable.el2z:
            raise ValidationError(
                """Position symbol not in periodic table: {}""".format(
                    self.__el))
        if self.__iso:
            if self.__el + self.__iso not in periodictable.eliso2masses:
                raise ValidationError(
                    """Position isotope not in periodic """
                    """table: {}""".format(
                        el + iso))

    def __eliso(self):
        """element label with isotope appropriate for periodic table key
        (e.g., CO, CO59, CO60)

        """
        if self.__iso:
            return self.__el + self.__iso
        else:
            return self.__el

    @property
    def Z(self):
        """nuclear charge for position"""
        return periodictable.el2z[self.__el]

    @property
    def mass(self):
        """mass for position"""
        try:
            return self.__mass
        except AttributeError:
            return periodictable.eliso2masses[self.__eliso()]

    @mass.setter
    def mass(self, val):
        if val is not None:
            self.__mass = float(val)

    @property
    def symbol(self):
        """element label for position, cleaned-up (C2 => C, CO_iii = Co)"""
        return self.__el.title()

    @property
    def charge(self):
        """Nuclear charge on location in electron units"""
        try:
            return self.__charge
        except AttributeError:
            return float(self.Z)

    @charge.setter
    def charge(self, val):
        if val is not None:
            self.__charge = float(val)


class MoleculeGist(object):
    """

    """

    def __init__(self, atoms):
        """

        """
        #: atom info vector
        self.atoms = []


if __name__ == '__main__':
    import sys
    sys.path.append('/Users/loriab/linux/qcdb/qcdb')

    print AtomGist('PU_123', 1.0, 2, -5)
    print AtomGist('c', 0.0, 0, 0, charge=-6)
#    a1 = AtomGist('pu_123', 0, 2, -5, charge=+95)
#    print a1
#    a1.x = 6 #'abc'
#    #a1.symbol = 'co'
#    a1.label = 'co'
#    print a1
    print AtomGist('pu123', 0.0, 0.0, 0.0)
    print AtomGist('239pu_abc', 0.0, 0.0, 0.0)
    a1 = AtomGist('239pu', 0.0, 0.0, 0.0)
    print dir()
    print a1.mass
