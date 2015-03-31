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
        #return """  %-23s %16.8f %16.8f %16.8f    %6.2f    %8.3f %s""" % (
        #    self.label + ' / ' + self.symbol + ' / ' + str(self.Z),
        #    self.x, self.y, self.z, self.charge, self.mass,
        #    ', '.join(['{}={}'.format(k, v) for k, v in self.tags.iteritems()]))
        return """  %-23s %16.8f %16.8f %16.8f    %8.3f""" % (
            self.label + ' / ' + self.symbol + ' / ' + str(self.Z),
            self.x, self.y, self.z, self.mass) #, self.charge, self.mass,
            #', '.join(['{}={}'.format(k, v) for k, v in self.tags.iteritems()]))

    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result

    @property
    def label(self):
        """label for position"""
        #ans = self.__iso if self.__iso else ''
        #ans += self.symbol
        #ans += self.__lbl if self.__lbl else ''
        ans = ''
        #try:
        ##if hasattr(self, '__iso'):
        #    ans += str(self.__iso)
        #    print('QWER', ans, self.__iso)
        #except AttributeError:
        #    pass
        ans += '' if self.isotope is None else str(self.isotope)
        ans += self.symbol
#        ans += self.__lbl if self.__lbl else ''
        #if hasattr(self, '__lbl'):
        try:
            ans += self.__lbl
            #print('ASDF', ans, self.__lbl)
        except AttributeError:
            pass
        return ans

    @label.setter
    def label(self, val):
        atm = atom.match(val.strip().capitalize())
        iso = atm.group('isotope')
        el = atm.group('symbol').upper()
        lbl = atm.group('label')
        #self.__iso = atm.group('isotope')
        #self.__el = atm.group('symbol').upper()
        print('QW', lbl)
        if lbl:
            self.__lbl = lbl
            print('QW2', self.__lbl)
        else:
            try:
                print('deleting __lbl', self.__lbl)
                del self.__lbl
            except AttributeError:
                pass
#        else:
#            self.__lbl = None
        #print('==>', self.__iso, self.__el, self.__lbl, '<==')
        print('==>', iso, el, lbl, '<==')
        #print(periodictable.eliso2masses['CL'])
        #print(periodictable.eliso2masses['CL35'])
        #print(periodictable.eliso2masses['CL37'])
        #print(periodictable.el2masses['CL'])
        if not atm:
            raise ValidationError(
                """Position label not of Co, Co3, 60Co (isotope), or """
                """Co_III form, case insensitive: {}.""".format(
                    val))

        if el in periodictable.el2z:
            self.Z = periodictable.el2z[el]
        else:
            raise ValidationError(
                """Position symbol not in periodic table: {}""".format(
                    el))

        if iso:
            if el + iso in periodictable.eliso2masses:
                self.__iso = iso
                print('setting iso', self.__iso)
            else:
                raise ValidationError(
                    """Position isotope not in periodic """
                    """table: {}""".format(
                       el + iso))
        try:
            print('NM', self.__lbl)
        except AttributeError:
            pass

    def __eliso(self):
        """element label with isotope appropriate for periodic table key
        (e.g., CO, CO59, CO60)

        """
        try:
            print('in eliso trying', self.symbol.upper(), str(self.__iso))
            return self.symbol.upper() + str(self.__iso)
        except AttributeError:
            return self.symbol.upper()
        #if self.__iso:
        #    return self.symbol.upper() + self.__iso
        #else:
        #    return self.symbol.upper()

    @property
    def Z(self):
        """nuclear charge for position"""
        #return periodictable.el2z[self.__el]
        return self.__z

    @Z.setter
    def Z(self, val):
        try:
            oldZ = self.__z
        except AttributeError:
            oldZ = -1
        #oldZ = self.__z if hasattr(self, '__z') else -1
        if val in periodictable.z2el.keys():
            self.__z = val
            if val != oldZ:
#                self.__mass = None
                self.__charge = None
                #self.__iso = None
                try:
                    del self.__iso
                    del self.__mass
                    del self.__lbl
                except AttributeError:
                    pass
                print('transmutation!', val, oldZ)
            else:
                print('no transmutation')
    
    @property
    def isotope(self):
        try:
            return self.__iso
        except AttributeError:
            return int(round(self.mass, 0))

    @property
    def symbol(self):
        """element label for position, cleaned-up (C2 => C, CO_iii = Co)"""
        return periodictable.z2el[self.__z].title()
            
    @property
    def mass(self):
        """mass for position"""
        try:
            return self.__mass
        except AttributeError:
            print('setting mass based on', self.__eliso())
            return periodictable.eliso2masses[self.__eliso()]

    @mass.setter
    def mass(self, val):
        if val is not None:
            fval = float(val)
            self.__mass = fval
            iso = int(round(fval, 0))
            eliso = self.symbol.upper() + str(iso)
            print(eliso, fval, periodictable.eliso2masses[eliso] - fval, 'TYTY')
            if abs(periodictable.eliso2masses[eliso] - fval) < 1.0e-3:
                self.__iso = iso
                print('setting isotope', iso)
            else:
                self.__iso = None
                #del self.__iso
                print('un-setting isotope', iso)

#    @property
#    def symbol(self):
#        """element label for position, cleaned-up (C2 => C, CO_iii = Co)"""
#        return self.__el.title()

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
    from psiutil import *

    print(AtomGist('PU_123', 1.0, 2, -5))
    print(AtomGist('c', 0.0, 0, 0, charge=-6))
#    a1 = AtomGist('pu_123', 0, 2, -5, charge=+95)
#    print a1
#    a1.x = 6 #'abc'
#    #a1.symbol = 'co'
#    a1.label = 'co'
#    print a1
    print(AtomGist('pu123', 0.0, 0.0, 0.0))
    print(AtomGist('239pu_abc', 0.0, 0.0, 0.0))
    a1 = AtomGist('239pu', 0.0, 0.0, 0.0)
    print(a1.mass)
    print(AtomGist('60co', 0.0, 0.0, 0.0))
    print(AtomGist('59co', 0.0, 0.0, 0.0))

    def checkAtom(expect, obj):
        if 'mass' in expect:
            attr = 'mass'
            compare_values(expect[attr], getattr(obj, attr), 3, '{} is {}'.format(attr, getattr(obj, attr)))
        if 'Z' in expect:
            attr = 'Z'
            compare_integers(expect[attr], getattr(obj, attr), '{} is {}'.format(attr, getattr(obj, attr)))
        if 'label' in expect:
            attr = 'label'
            compare_strings(expect[attr], getattr(obj, attr), '{} is {}'.format(attr, getattr(obj, attr)))
    print('      Initial plain atom')



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
    import copy
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
    #expect = {'mass': , 'Z': , 'label': }
    #checkAtom(expect, a2)


#    print("""      Set isotope via label: label = '60co_aBc'""")
#    a2.label = '60co_aBc'
#    expect = {'mass': 59.933817059, 'Z': 27, 'label': '60Co_abc'}
#    checkAtom(expect, a2)
#
#    print('      Set isotope via mass: mass = 61')
#    a2.mass = 61
#    expect = {'mass': 61, 'Z': 27, 'label': 'Co_abc'}
#    checkAtom(expect, a2)
#
#    print("""      Set label: label = '59CO_Abc'""")
#    a2.label = '59Co_Abc'
#    print(a2.mass)
#    print(a2.Z)
#    print(a2.label)



    #attr = 'mass'; ans = 61; compare_values(ans, getattr(a2, attr), 5, '{} is {}'.format(attr, ans))
    #attr = 'Z'; ans = 27; compare_integers(ans, getattr(a2, attr), '{} is {}'.format(attr, ans))
    #attr = 'label'; ans = 'Co_abc'; compare_strings(ans, getattr(a2, attr), '{} is {}'.format(attr, ans))

