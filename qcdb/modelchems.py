

class QCEssential(object):
    """Class to link literature and external representation of some
    aspect of quantum chemistry (basis set, method, etc.) with a
    shorthand and indexed representation of same.

    """
    def __init__(self, name, fullname=None, latex=None, citation=None, pdfdatabase=None, comment=None):
        """

        """
        self.name = name.lower()
        self.fullname = fullname
        if fullname is not None and latex is None:
            self.latex = fullname
        else:
            self.latex = latex
        self.citation = citation
        self.dsdbid = pdfdatabase
        self.comment = comment

    def __str__(self):
        text = ''
        text += """  ==> %s BasisSet <==\n\n""" % (self.name)
        text += """  Formal name:          %s\n""" % (self.fullname)
        text += """  LaTeX representation: %s\n""" % (self.latex)
        text += """  PDF database id:      %s\n""" % (self.dsdbid)
        text += """  Literature citation:  %s\n""" % (self.citation)
        text += """  Comment:              %s\n""" % (self.comment)
        text += """\n"""
        return text


class BasisSet(QCEssential):
    """Specialization of :pyclass:`QCEssential` for basis sets.

    """
    def __init__(self, name, fullname=None, latex=None, citation=None, pdfdatabase=None, comment=None, zeta=None):
        QCEssential.__init__(self, name, fullname, latex, citation, pdfdatabase, comment)
        self.name = name.lower()
        self.zeta = zeta

    def __str__(self):
        text = ''
        text += """  ==> %s BasisSet <==\n\n""" % (self.name)
        text += """  Formal name:          %s\n""" % (self.fullname)
        text += """  LaTeX representation: %s\n""" % (self.latex)
        text += """  Zeta:                 %s\n""" % (self.zeta)
        text += """  PDF database id:      %s\n""" % (self.dsdbid)
        text += """  Literature citation:  %s\n""" % (self.citation)
        text += """  Comment:              %s\n""" % (self.comment)
        text += """\n"""
        return text


class Method(QCEssential):
    """Specialization of :pyclass:`QCEssential` for quantum chemical methods.

    """
    def __init__(self, name, fullname=None, latex=None, citation=None, pdfdatabase=None, comment=None):
        QCEssential.__init__(self, name, fullname, latex, citation, pdfdatabase, comment)
        self.name = name.upper()

    def __str__(self):
        text = ''
        text += """  ==> %s Method <==\n\n""" % (self.name)
        text += """  Formal name:          %s\n""" % (self.fullname)
        text += """  LaTeX representation: %s\n""" % (self.latex)
        text += """  PDF database id:      %s\n""" % (self.dsdbid)
        text += """  Literature citation:  %s\n""" % (self.citation)
        text += """  Comment:              %s\n""" % (self.comment)
        text += """\n"""
        return text


_tlist = [
    BasisSet('dz',         fullname='cc-pVDZ'),
    BasisSet('jadz',       fullname='jun-cc-pVDZ'),
    BasisSet('hadz',       fullname='heavy-aug-cc-pVDZ'),
    BasisSet('adz',        fullname='aug-cc-pVDZ'),
    BasisSet('tz',         fullname='cc-pVTZ'),
    BasisSet('matz',       fullname='may-cc-pVTZ'),
    BasisSet('jatz',       fullname='jun-cc-pVTZ'),
    BasisSet('hatz',       fullname='heavy-aug-cc-pVTZ'),
    BasisSet('atz',        fullname='aug-cc-pVTZ'),
    BasisSet('qz',         fullname='cc-pVQZ'),
    BasisSet('aaqz',       fullname='apr-cc-pVQZ'),
    BasisSet('maqz',       fullname='may-cc-pVQZ'),
    BasisSet('jaqz',       fullname='jun-cc-pVQZ'),
    BasisSet('haqz',       fullname='heavy-aug-cc-pVQZ'),
    BasisSet('aqz',        fullname='aug-cc-pVQZ'),
    BasisSet('a5z',        fullname='aug-cc-pV5Z'),
    BasisSet('dtz',        fullname='cc-pVDTZ'),
    BasisSet('jadtz',      fullname='jun-cc-pVDTZ'),
    BasisSet('hadtz',      fullname='heavy-aug-cc-pVDTZ'),
    BasisSet('adtz',       fullname='aug-cc-pVDTZ'),
    BasisSet('tqz',        fullname='cc-pVTQZ'),
    BasisSet('matqz',      fullname='may-cc-pVTQZ'),
    BasisSet('jatqz',      fullname='jun-cc-pVTQZ'),
    BasisSet('hatqz',      fullname='heavy-aug-cc-pVTQZ'),
    BasisSet('atqz',       fullname='aug-cc-pVTQZ'),
    BasisSet('aq5z',       fullname='aug-cc-pVQ5Z'),
    BasisSet('a6z',        fullname='aug-cc-pV6Z'),
    BasisSet('a56z',       fullname='aug-cc-pV56Z'),
    BasisSet('atzdz',      fullname='[aTZ; D:DZ]', latex="""[aTZ; $\delta$:DZ]"""),
    BasisSet('adtzdz',     fullname='[aDTZ; D:DZ]', latex="""[aDTZ; $\delta$:DZ]"""),
    BasisSet('atqzdz',     fullname='[aTQZ; D:DZ]', latex="""[aTQZ; $\delta$:DZ]"""),
    BasisSet('atzjadz',    fullname='[aTZ; D:jaDZ]', latex="""[aTZ; $\delta$:jaDZ]"""),
    BasisSet('adtzjadz',   fullname='[aDTZ; D:jaDZ]', latex="""[aDTZ; $\delta$:jaDZ]"""),
    BasisSet('atqzjadz',   fullname='[aTQZ; D:jaDZ]', latex="""[aTQZ; $\delta$:jaDZ]"""),
    BasisSet('atzhadz',    fullname='[aTZ; D:haDZ]', latex="""[aTZ; $\delta$:haDZ]"""),
    BasisSet('adtzhadz',   fullname='[aDTZ; D:haDZ]', latex="""[aDTZ; $\delta$:haDZ]"""),
    BasisSet('atqzhadz',   fullname='[aTQZ; D:haDZ]', latex="""[aTQZ; $\delta$:haDZ]"""),
    BasisSet('atzadz',     fullname='[aTZ; D:aDZ]', latex="""[aTZ; $\delta$:aDZ]"""),
    BasisSet('adtzadz',    fullname='[aDTZ; D:aDZ]', latex="""[aDTZ; $\delta$:aDZ]"""),
    BasisSet('atqzadz',    fullname='[aTQZ; D:aDZ]', latex="""[aTQZ; $\delta$:aDZ]"""),
    BasisSet('aq5zadz',    fullname='[aQ5Z; D:aDZ]', latex="""[aQ5Z; $\delta$:aDZ]"""),
    BasisSet('atzdtz',     fullname='[aTZ; D:DTZ]', latex="""[aTZ; $\delta$:DTZ]"""),
    BasisSet('atqzdtz',    fullname='[aTQZ; D:DTZ]', latex="""[aTQZ; $\delta$:DTZ]"""),
    BasisSet('atzjadtz',   fullname='[aTZ; D:jaDTZ]', latex="""[aTZ; $\delta$:jaDTZ]"""),
    BasisSet('atqzjadtz',  fullname='[aTQZ; D:jaDTZ]', latex="""[aTQZ; $\delta$:jaDTZ]"""),
    BasisSet('atzhadtz',   fullname='[aTZ; D:haDTZ]', latex="""[aTZ; $\delta$:haDTZ]"""),
    BasisSet('atqzhadtz',  fullname='[aTQZ; D:haDTZ]', latex="""[aTQZ; $\delta$:haDTZ]"""),
    BasisSet('atzadtz',    fullname='[aTZ; D:aDTZ]', latex="""[aTZ; $\delta$:aDTZ]"""),
    BasisSet('atqzadtz',   fullname='[aTQZ; D:aDTZ]', latex="""[aTQZ; $\delta$:aDTZ]"""),
    BasisSet('aq5zadtz',   fullname='[aQ5Z; D:aDTZ]', latex="""[aQ5Z; $\delta$:aDTZ]"""),
    BasisSet('atqztz',     fullname='[aTQZ; D:TZ]', latex="""[aTQZ; $\delta$:TZ]"""),
    BasisSet('atqzjatz',   fullname='[aTQZ; D:jaTZ]', latex="""[aTQZ; $\delta$:jaTZ]"""),
    BasisSet('atqzmatz',   fullname='[aTQZ; D:maTZ]', latex="""[aTQZ; $\delta$:maTZ]"""),
    BasisSet('atqzhatz',   fullname='[aTQZ; D:haTZ]', latex="""[aTQZ; $\delta$:haTZ]"""),
    BasisSet('atqzatz',    fullname='[aTQZ; D:aTZ]', latex="""[aTQZ; $\delta$:aTZ]"""),
    BasisSet('aq5zatz',    fullname='[aQ5Z; D:aTZ]', latex="""[aQ5Z; $\delta$:aTZ]"""),
    ]
bases = {}
for item in _tlist:
    bases[item.name] = item

_tlist = [
    Method('SAPT0',           fullname='SAPT0'),
    Method('SAPT0S',          fullname='sSAPT0', latex="""\\textit{s}SAPT0"""),
    Method('SAPTDFT',         fullname='DFT-SAPT'),
    Method('SAPT2',           fullname='SAPT2'),
    Method('SAPT2P',          fullname='SAPT2+'),
    Method('SAPT3',           fullname='SAPT2+(3)'),
    Method('SAPT3F',          fullname='SAPT2+3'),
    Method('SAPT2PC',         fullname='SAPT2+(CCD)'),
    Method('SAPT3C',          fullname='SAPT2+(3)(CCD)'),
    Method('SAPTCCD',         fullname='SAPT2+3(CCD)'),
    Method('SAPT2PM',         fullname='SAPT2+dMP2', latex="""SAPT2+$\delta$MP2"""),
    Method('SAPT3M',          fullname='SAPT2+(3)dMP2', latex="""SAPT2+(3)$\delta$MP2"""),
    Method('SAPT3FM',         fullname='SAPT2+3dMP2', latex="""SAPT2+3$\delta$MP2"""),
    Method('SAPT2PCM',        fullname='SAPT2+(CCD)dMP2', latex="""SAPT2+(CCD)$\delta$MP2"""),
    Method('SAPT3CM',         fullname='SAPT2+(3)(CCD)dMP2', latex="""SAPT2+(3)(CCD)$\delta$MP2"""),
    Method('SAPTCCDM',        fullname='SAPT2+3(CCD)dMP2', latex="""SAPT2+3(CCD)$\delta$MP2"""),
    Method('SAPT2LCM',        fullname='MP2(CCD)', comment="""Identical to SAPT2+(CCD)dMP2"""),
    Method('HF',              fullname='HF'),
    Method('MP2',             fullname='MP2'),
    Method('SCSMP2',          fullname='SCS-MP2'),
    Method('SCSNMP2',         fullname='SCS(N)-MP2'),
    Method('SCSMIMP2',        fullname='SCS(MI)-MP2'),
    Method('DWMP2',           fullname='DW-MP2'),
    Method('MP2C',            fullname='MP2C'),
    Method('MP3',             fullname='MP3'),
    Method('MP25',            fullname='MP2.5'),
    Method('CCSD',            fullname='CCSD'),
    Method('SCSCCSD',         fullname='SCS-CCSD'),
    Method('SCSMICCSD',       fullname='SCS(MI)-CCSD'),
    Method('CCSDT',           fullname='CCSD(T)'),
    Method('HFCABS',          fullname='HF-CABS'),
    Method('MP2F12',          fullname='MP2-F12'),
    Method('SCSMP2F12',       fullname='SCS-MP2-F12'),
    Method('SCSNMP2F12',      fullname='SCS(N)-MP2-F12'),
    Method('SCSMIMP2F12',     fullname='SCS(MI)-MP2-F12'),
    Method('DWMP2F12',        fullname='DW-MP2-F12'),
    Method('MP2CF12',         fullname='MP2C-F12'),
    Method('CCSDAF12',        fullname='CCSD-F12a'),
    Method('CCSDBF12',        fullname='CCSD-F12b'),
    Method('SCSCCSDAF12',     fullname='SCS-CCSD-F12a'),
    Method('SCSCCSDBF12',     fullname='SCS-CCSD-F12b'),
    Method('SCMICCSDAF12',    fullname='SCS(MI)-CCSD-F12a'),
    Method('SCMICCSDBF12',    fullname='SCS(MI)-CCSD-F12b'),
    Method('CCSDTAF12',       fullname='CCSD(T)-F12a'),
    Method('CCSDTBF12',       fullname='CCSD(T)-F12b'),
    Method('DWCCSDTF12',      fullname='DW-CCSD(T)-F12'),
    Method('B97D3',           fullname='B97-D3'),
    Method('B3LYPD3',         fullname='B3LYP-D3'),
    Method('wB97XD',          fullname='$\omega$B97X-D'),
    Method('M062X',           fullname='M06-2X'),
    ]

methods = {}
for item in _tlist:
    methods[item.name] = item
