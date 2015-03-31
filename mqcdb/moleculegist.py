"""

"""
from __future__ import print_function
import re
import periodictable
from xcpt import *

atom = re.compile(r"""^(?P<isotope>[0-9]{1,3})?  # optional isotope
                       (?P<symbol>[A-Z]{1,3})    # mandatory atomic symbol
                       (?P<label>(_\w+)|(\d+))?$ # optional userid""",
                  re.IGNORECASE | re.X)


class AtomGist(object):
    """

    """

    def __init__(self, label, x, y, z, mass=None, charge=None, tags={}):
        """

        """
        #: atomic number
        self.__z = None

        #: mass in Daltons
        self.__mass = None

        #:
        self.__charge = None

        #: user part of label
        self.__lbl = ''

        #: x in Bohr
        self.x = float(x)

        #: y in Bohr
        self.y = float(y)

        #: z in Bohr
        self.z = float(z)

        # constructor
        self.label = label
        self.mass = mass
        self.charge = charge

        #: all tags associated with position
        self.tags = {}
        for item, val in tags.iteritems():
            self.tags[item.lower()] = val

    def __str__(self):
        #return """  %-23s %16.8f %16.8f %16.8f    %6.2f    %8.3f %s""" % (
        #    self.label + ' / ' + self.symbol + ' / ' + str(self.Z),
        #    self.x, self.y, self.z, self.charge, self.mass,
        #    ', '.join(['{}={}'.format(k, v) for k, v in self.tags.iteritems()]))
        return """  %-23s %16.8f %16.8f %16.8f    %8.3f""" % (
            self.label + ' / ' + self.symbol + ' / ' + str(self.Z),
            self.x, self.y, self.z, self.mass)

#    def __copy__(self):
#        cls = self.__class__
#        result = cls.__new__(cls)
#        result.__dict__.update(self.__dict__)
#        return result

    @property
    def label(self):
        """label for position"""
        return self.nuclide() + self.__lbl

    @label.setter
    def label(self, val):
        atm = atom.match(val.strip().lower())
        if not atm:
            raise ValidationError(
                """Position label not of Co, Co3, 60Co (nuclide), or """
                """Co_III form, case insensitive: {}.""".format(
                    val))
        self.symbol = atm.group('symbol')

        iso = atm.group('isotope')
        if iso:
            self.A = iso

        lbl = atm.group('label')
        if lbl:
            self.__lbl = lbl.lower()
        else:
            self.__lbl = ''

    @property
    def Z(self):
        """nuclear charge for position"""
        return self.__z

    @Z.setter
    def Z(self, val):
        oldZ = self.__z
        val = int(val)
        if val in periodictable.z2el.keys():
            self.__z = val
            if val != oldZ:
                self.__mass = None
                self.__charge = None
                self.__lbl = ''

    @property
    def symbol(self):
        """element label for position, cleaned-up (C2 => C, CO_iii = Co)"""
        return periodictable.z2el[self.__z].title()

    @symbol.setter
    def symbol(self, val):
        val = val.upper()
        if val in periodictable.el2z:
            self.Z = periodictable.el2z[val]
        else:
            raise ValidationError(
                """Position symbol not in periodic table: {}""".format(
                    val))

    @property
    def mass(self):
        """mass for position"""
        if self.__mass:
            return self.__mass
        else:
            return periodictable.z2mass[self.Z]

    @mass.setter
    def mass(self, val):
        if val is not None:
            self.__mass = float(val)

    @property
    def A(self):
        """nucleon number for position"""
        nucleon = int(round(self.mass, 0))
        libMass = periodictable.eliso2mass[self.symbol.upper() + str(nucleon)]
        if abs(libMass - self.mass) < 1.0e-3:
            return nucleon

    @A.setter
    def A(self, val):
        eliso = self.symbol.upper() + str(val)
        if eliso in periodictable.eliso2mass:
            self.mass = periodictable.eliso2mass[eliso]
        else:
            raise ValidationError(
                """Position isotope not in periodic table: {}""".format(
                    el + iso))

    def nuclide(self):
        """element label with isotope appropriate for periodic table key
        (e.g., Co, 59Co, 60Co)

        """
        if self.A:
            return str(self.A) + self.symbol
        else:
            return self.symbol

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
    import copy
    sys.path.append('/Users/loriab/linux/qcdb/qcdb')
    from psiutil import *

    def checkAtom(expect, obj):
        if 'mass' in expect:
            attr = 'mass'
            compare_values(expect[attr], getattr(obj, attr), 3,
                           '{} is {}'.format(attr, getattr(obj, attr)))
        if 'Z' in expect:
            attr = 'Z'
            compare_integers(expect[attr], getattr(obj, attr),
                             '{} is {}'.format(attr, getattr(obj, attr)))
        if 'label' in expect:
            attr = 'label'
            compare_strings(expect[attr], getattr(obj, attr),
                            '{} is {}'.format(attr, getattr(obj, attr)))

    print("""\n      AtomGist('co_abC', 0.0, 0.0, 0.0)""")
    a2 = AtomGist('co_abC', 0.0, 0.0, 0.0)
    expect = {'mass': 58.933195048, 'Z': 27, 'label': '59Co_abc'}
    checkAtom(expect, a2)

    print("""\n      a2.label = '60co_abC'""")
    a2.label = '60co_abC'
    expect = {'mass': 59.933817059, 'Z': 27, 'label': '60Co_abc'}
    checkAtom(expect, a2)

    print("""\n      a2.mass = 61""")
    a2.mass = 61
    expect = {'mass': 61.0, 'Z': 27, 'label': 'Co_abc'}
    checkAtom(expect, a2)

    print("""\n      a3 = copy.copy(a2)\n      a2.label = 'u'""")
    a3 = copy.copy(a2)
    # resetting Z wipes out all mass, charge info
    a2.label = 'u'
    expect = {'mass': 238.050788247, 'Z': 92, 'label': '238U'}
    checkAtom(expect, a2)

    print("""\n      check a3""")
    # check copy is independent
    expect = {'mass': 61.0, 'Z': 27, 'label': 'Co_abc'}
    checkAtom(expect, a3)

    print("""\n      a3.Z = 92""")
    a3.Z = 92
    expect = {'mass': 238.050788247, 'Z': 92, 'label': '238U'}
    checkAtom(expect, a3)

    print("""\n      a3.mass = 235.044""")
    a3.mass = 235.044
    expect = {'mass': 235.043929918, 'Z': 92, 'label': '235U'}
    checkAtom(expect, a3)

    print("""\n      a3.label = '238U'""")
    a3.label = '238U'
    expect = {'mass': 238.050788247, 'Z': 92, 'label': '238U'}
    checkAtom(expect, a3)

    print("""\n      a3.label = '234U_trace'""")
    a3.label = '234U_trace'
    expect = {'mass': 234.040952088, 'Z': 92, 'label': '234U_trace'}
    checkAtom(expect, a3)

    print("""\n      a3.mass = 234""")
    a3.mass = 234
    expect = {'mass': 234.0, 'Z': 92, 'label': 'U_trace'}
    checkAtom(expect, a3)

    # expect = {'mass': , 'Z': , 'label': }
    # checkAtom(expect, a2)
