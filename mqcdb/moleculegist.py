"""

"""
from __future__ import print_function
import re
import itertools
import periodictable
from xcpt import *

atom = re.compile(r"""^(?P<isotope>[0-9]{1,3})?  # optional isotope
                       (?P<symbol>[A-Z]{1,3})    # mandatory atomic symbol
                       (?P<label>(_\w+)|(\d+))?$ # optional userid""",
                  re.IGNORECASE | re.X)


class AtomGist(object):
    """Stores essential characteristics of a location, and provides basic
    interaction with them.

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
        label = '{0:>18s}'.format(
            self.label + ' / ' + self.symbol + ' / ' + str(self.Z))
        xyz = '{0:16.8f} {1:16.8f} {2:16.8f}'.format(
            self.x, self.y, self.z)
        mass = '{0:8.3f}'.format(self.mass)
        tags = '{0}'.format(', '.join('{}={}'.format(k, v) for k, v in
                            self.tags.iteritems()) if self.tags else '')
        return """  {0} {1} {2} {3}""".format(label, xyz, mass, tags)

        #return ("""  {l:>18s} {x:16.8f} {y:16.8f} {z:16.8f}"""
        #        """ {m:8.3f} {t}""".format(
        #           l=self.label + ' / ' + self.symbol + ' / ' + str(self.Z),
        #           x=self.x, y=self.y, z=self.z,
        #           m=self.mass,
        #           t=', '.join('{}={}'.format(k, v) for k, v in
        #                       self.tags.iteritems()) if self.tags else ''))

    #def __repr__(self):
    #    text = """AtomGist(
    #def __init__(self, label, x, y, z, mass=None, charge=None, tags={}):


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

    def __init__(self, atoms, name=None):
        """

        """
        #: Molecule name
        self.__name = None

        #: atom info vector
        self.__atoms = []

        self.name = name
        self.atoms = atoms

    @property
    def name(self):
        """label for object"""
        return '' if self.__name is None else self.__name

    @name.setter
    def name(self, val):
        if val is not None:
            if val.isalnum():
                self.__name = val.lower()
            else:
                raise ValidationError(
                    """Molecule name not alphanumeric: {}""".format(
                        val))

    @property
    def atoms(self):
        """array of locations which Molecule comprises"""
        return self.__atoms

    @atoms.setter
    def atoms(self, val):
        for at in val:
            if isinstance(at, AtomGist):
                self.__atoms.append(at)
            else:
                raise ValidationError(
                    """Building Molecule not from AtomGist: {}""".format(
                        at))

    def __str__(self):
        return self.print_out()

    def print_out(self, angstrom=False, sset=None):
        atomGenerator = self.dissect_sset(sset)
        text = ''
        text += """  ==> {} MoleculeGist <==\n\n""".format(self.name)
        text += """  Geometry (in {}):\n""".format(
            'Angstrom' if angstrom else 'Bohr')
        for at in atomGenerator:
            text += at.__str__() + '\n'
        return text

    def dissect_sset(self, sset):
        """"""
        if sset in ['all', 'real']:
            if sset == 'all':
                return (at for at in self.atoms)
            elif sset == 'real':
                return (at for at in self.atoms if at.Z > 0)
                # TODO not right test
        else:
            return itertools.ifilter(sset, self.atoms)

    def generateAllAtom(self):
        for at in self.atoms:
            yield at

    def hastag(self, tag):
        """Returns generator of atoms with *tag*"""
        for at in self.atoms:
            if tag.lower() in at.tags:
                yield at

#    def real(self):
#        """"""
#        for at in self.atoms:
#            if at.Z > 0:  # TODO not right test
#                yield at

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

    ams = [
        AtomGist('O', 0.0, 0.0, 0.0),
        AtomGist('H', 1.0, 0.0, 0.0, tags={'star': True, 'galaxy': 'five'}),
        AtomGist('H', 0.0, 1.0, 0.0, tags={'star': False}),
        AtomGist('X', 1.0, 0.0, 3.0),
        AtomGist('H_dim', 0.0, 1.0, 3.0),
        AtomGist('H_dim', 0.0, 1.0, 3.0)]

    mol1 = MoleculeGist(ams)
    mol1.name = 'mymol'
    print(mol1.print_out())
    print(mol1.print_out(sset='real'))
    print(mol1.print_out(sset=lambda a: a.label.endswith('H_dim') or a.Z == 8))
    print(mol1.print_out(sset=lambda a: 'star' in a.tags))
    print(mol1.print_out(
        sset=lambda a: 'star' in a.tags and a.tags['star'] is False))
