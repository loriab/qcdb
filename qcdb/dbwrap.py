import os
import sys
import math
import itertools
import collections
from exceptions import *
from modelchems import Method, BasisSet, Error, methods, bases, errors
import psiutil
import textables


def initialize_errors(e=None, pe=None, pbe=None, extrema=True):
    """

    """
    error = collections.OrderedDict()
    error['maxe'] = None if (e is None or not extrema) else e        # LD_XA
    error['mine'] = None if (e is None or not extrema) else e        # LD_XI
    error['me'] = None if e is None else 0.0                         # LD_MS
    error['mae'] = None if e is None else 0.0                        # LD_MA
    error['rmse'] = None if e is None else 0.0                       # LD_RA
    error['stde'] = None if e is None else 0.0
    error['maxpe'] = None if (pe is None or not extrema) else pe     # FD_XA
    error['minpe'] = None if (pe is None or not extrema) else pe     # FD_XI
    error['mpe'] = None if pe is None else 0.0                       # FD_MS
    error['mape'] = None if pe is None else 0.0                      # FD_MA
    error['rmspe'] = None if pe is None else 0.0                     # FD_RA
    error['stdpe'] = None if pe is None else 0.0
    error['maxpbe'] = None if (pbe is None or not extrema) else pbe  # BD_XA
    error['minpbe'] = None if (pbe is None or not extrema) else pbe  # BD_XI
    error['mpbe'] = None if pbe is None else 0.0                     # BD_MS
    error['mapbe'] = None if pbe is None else 0.0                    # BD_MA
    error['rmspbe'] = None if pbe is None else 0.0                   # BD_RA
    error['stdpbe'] = None if pbe is None else 0.0
    return error


def average_errors(*args):
    """Each item in *args* should be an error dictionary. Performs
    average-like operation over all items, which should be error
    dictionaries, in *args*. Defined for ME, MAE, STDE, and their
    relative-error variants. None returned for undefined statistics or
    when an item is missing.

    """
    Ndb = float(len(args))
    avgerror = initialize_errors()
    try:
        avgerror['me'] = sum([x['me'] for x in args]) / Ndb
        avgerror['mae'] = sum([x['mae'] for x in args]) / Ndb
        avgerror['stde'] = math.sqrt(sum([x['stde'] * x['stde'] for x in args]) / Ndb)
        avgerror['mpe'] = sum([x['mpe'] for x in args]) / Ndb
        avgerror['mape'] = sum([x['mape'] for x in args]) / Ndb
        avgerror['stdpe'] = math.sqrt(sum([x['stdpe'] * x['stdpe'] for x in args]) / Ndb)
        avgerror['mpbe'] = sum([x['mpbe'] for x in args]) / Ndb
        avgerror['mapbe'] = sum([x['mapbe'] for x in args]) / Ndb
        avgerror['stdpbe'] = math.sqrt(sum([x['stdpbe'] * x['stdpbe'] for x in args]) / Ndb)
    except TypeError:
        pass
    return avgerror


def format_errors(err, mode=1):
    """From error dictionary *err*, returns a LaTeX-formatted string,
    after handling None entries.

    """
    if mode == 1:
        me = ' ----' if err['me'] is None else '%+.2f' % (err['me'])
        stde = '----' if err['stde'] is None else '%.2f' % (err['stde'])
        mae = '  ----' if err['mae'] is None else '%6.2f' % (err['mae'])
        mape = '  ----  ' if err['mape'] is None else '%6.1f\%%' % (100 * err['mape'])
        mapbe = '  ----  ' if err['mapbe'] is None else '%6.1f\%%' % (100 * err['mapbe'])
        text = """$[%s, %s]$ %s %s %s""" % \
            (me, stde, mae, mape, mapbe)
        return text

    if mode == 2:
        maxe = '----' if err['maxe'] is None else '%8.2f' % (err['maxe'])
        mine = '----' if err['mine'] is None else '%8.2f' % (err['mine'])
        me = '----' if err['me'] is None else '%+8.2f' % (err['me'])
        mae = '----' if err['mae'] is None else '%8.2f' % (err['mae'])
        rmse = '----' if err['rmse'] is None else '%8.2f' % (err['rmse'])
        stde = '----' if err['stde'] is None else '%8.2f' % (err['stde'])
        maxpe = '----' if err['maxpe'] is None else '%8.1f' % (100 * err['maxpe'])
        minpe = '----' if err['minpe'] is None else '%8.1f' % (100 * err['minpe'])
        mpe = '----' if err['mpe'] is None else '%+8.1f' % (100 * err['mpe'])
        mape = '----' if err['mape'] is None else '%8.1f' % (100 * err['mape'])
        rmspe = '----' if err['rmspe'] is None else '%8.1f' % (100 * err['rmspe'])
        stdpe = '----' if err['stdpe'] is None else '%8.1f' % (100 * err['stdpe'])
        maxpbe = '----' if err['maxpbe'] is None else '%8.1f' % (100 * err['maxpbe'])
        minpbe = '----' if err['minpbe'] is None else '%8.1f' % (100 * err['minpbe'])
        mpbe = '----' if err['mpbe'] is None else '%+8.1f' % (100 * err['mpbe'])
        mapbe = '----' if err['mapbe'] is None else '%8.1f' % (100 * err['mapbe'])
        rmspbe = '----' if err['rmspbe'] is None else '%8.1f' % (100 * err['rmspbe'])
        stdpbe = '----' if err['stdpbe'] is None else '%8.1f' % (100 * err['stdpbe'])
        text = """min: %s%s%s\nmax: %s%s%s\nm:   %s%s%s\nma:  %s%s%s\nrms: %s%s%s\nstd: %s%s%s""" % \
            (mine, minpe, minpbe, maxe, maxpe, maxpbe, me, mpe, mpbe, \
            mae, mape, mapbe, rmse, rmspe, rmspbe, stde, stdpe, stdpbe)
        return text

    if mode == 3:
        sdict = collections.OrderedDict()
        sdict['maxe'] = '' if err['maxe'] is None else '%8.2f' % (err['maxe'])
        sdict['mine'] = '' if err['mine'] is None else '%8.2f' % (err['mine'])
        sdict['me'] = '' if err['me'] is None else '%+8.2f' % (err['me'])
        sdict['mae'] = '' if err['mae'] is None else '%8.2f' % (err['mae'])
        sdict['rmse'] = '' if err['rmse'] is None else '%8.2f' % (err['rmse'])
        sdict['stde'] = '' if err['stde'] is None else '%8.2f' % (err['stde'])
        sdict['maxpe'] = '' if err['maxpe'] is None else '%8.1f' % (100 * err['maxpe'])
        sdict['minpe'] = '' if err['minpe'] is None else '%8.1f' % (100 * err['minpe'])
        sdict['mpe'] = '' if err['mpe'] is None else '%+8.1f' % (100 * err['mpe'])
        sdict['mape'] = '' if err['mape'] is None else '%8.1f' % (100 * err['mape'])
        sdict['rmspe'] = '' if err['rmspe'] is None else '%8.1f' % (100 * err['rmspe'])
        sdict['stdpe'] = '' if err['stdpe'] is None else '%8.1f' % (100 * err['stdpe'])
        sdict['maxpbe'] = '' if err['maxpbe'] is None else '%8.1f' % (100 * err['maxpbe'])
        sdict['minpbe'] = '' if err['minpbe'] is None else '%8.1f' % (100 * err['minpbe'])
        sdict['mpbe'] = '' if err['mpbe'] is None else '%+8.1f' % (100 * err['mpbe'])
        sdict['mapbe'] = '' if err['mapbe'] is None else '%8.1f' % (100 * err['mapbe'])
        sdict['rmspbe'] = '' if err['rmspbe'] is None else '%8.1f' % (100 * err['rmspbe'])
        sdict['stdpbe'] = '' if err['stdpbe'] is None else '%8.1f' % (100 * err['stdpbe'])
        return sdict


def string_contrast(ss):
    """From an array of strings, *ss*, returns maximum common prefix
    string, maximum common suffix string, and array of middles.

    """
    s = [item + 'q' for item in ss if item is not None]
    short = min(s, key=len)
    for ib in range(len(short)):
        if not all([mc[ib] == short[ib] for mc in s]):
            preidx = ib
            break
    else:
        preidx = 0
    for ib in range(len(short)):
        ie = -1 * (ib + 1)
        if not all([mc[ie] == short[ie] for mc in s]):
            sufidx = ie + 1
            break
    else:
        sufidx = -1 * (len(short))

    miditer = iter([mc[preidx:sufidx] for mc in s])
    prefix = short[:preidx]
    suffix = short[sufidx:-1]
    middle = ['' if mc is None else next(miditer) for mc in ss]

    return prefix, suffix, middle


#def expand_mc(modelchem_fragments, hierarchy):
#    """
#
#    """
#    modelchems = []
#    mc_hier = [modelchem_fragments[index] for index in hierarchy]
#    for combo in itertools.product(*mc_hier):
#        modelchems.append('-'.join([combo[hierarchy.index(ordinal)] for ordinal in range(len(hierarchy))]))
#    return modelchems


class ReactionDatum(object):
    """Piece of quantum chemical information that describes a qcdb.Reaction object.

    """
    def __init__(self, dbse, rxn, method, mode, basis, value, units='kcal/mol', comment=None):
        # geometry
        self.dbrxn = dbse + '-' + str(rxn)
        # qcdb.Method
        self.method = method
        # mode, e.g., unCP, CP, RLX, etc.
        self.mode = mode
        # qcdb.BasisSet
        self.basis = basis
        # numerical value for reaction
        self.value = float(value)
        # energy unit attached to value, defaults to kcal/mol
        self.units = units
        # addl comments
        self.comment = comment

    @classmethod
    def library_modelchem(cls, dbse, rxn, method, mode, basis, value, units='kcal/mol', comment=None):
        """Constructor when method and basis are strings corresponding to
        qcdb.Method and qcdb.BasisSet already defined in methods and bases.

        """
        # computational method
        if method.upper() in methods:
            tmp_method = methods[method.upper()]
        else:
            raise ValidationError("""Invalid ReactionDatum method %s.""" % (method))
        # computational basis set
        if basis.lower() in bases:
            tmp_basis = bases[basis.lower()]
        else:
            raise ValidationError("""Invalid ReactionDatum basis %s.""" % (basis))
        return cls(dbse, rxn, tmp_method, mode, tmp_basis, value, units='kcal/mol', comment=None)

    def __str__(self):
        text = ''
        text += """  ==> ReactionDatum <==\n\n"""
        text += """  Database reaction:    %s\n""" % (self.dbrxn)
        text += """  Method:               %s\n""" % (self.method.fullname)
        text += """  Mode:                 %s\n""" % (self.mode)
        text += """  Basis:                %s\n""" % (self.basis.fullname)
        text += """  Value:                %f [%s]\n""" % (self.value, self.units)
        text += """  Comment:              %s\n""" % (self.comment)
        text += """\n"""
        return text


class Reagent(object):
    """Chemical entity only slightly dresed up from qcdb.Molecule.

    """

    def __init__(self, name, mol, tagl=None, comment=None):
        # full name, e.g., 'S22-2-dimer' or 'NBC1-BzMe-8.0-monoA-CP' or 'HTBH-HCl-reagent'
        self.name = name
        # qcdb.Molecule
        try:
            self.NRE = mol.nuclear_repulsion_energy()
        except AttributeError:
            raise ValidationError("""Reagent must be instantiated with qcdb.Molecule object.""")
        else:
            self.mol = mol
        # description line
        self.tagl = tagl
        # addl comments
        self.comment = comment

    def __str__(self):
        text = ''
        text += """  ==> %s Reagent <==\n\n""" % (self.name)
        text += """  Tagline:              %s\n""" % (self.tagl)
        text += """  Comment:              %s\n""" % (self.comment)
        text += """  NRE:                  %f\n""" % (self.nre)
        text += """  Molecule:             %s\n"""
        text += """\n"""
        return text


class Reaction(object):
    """

    """

    def __init__(self, name, dbse, indx, tagl=None, latex=None, color='black', comment=None):
        # name, e.g., '2' or 'BzMe-8.0'
        self.name = name
        # database reaction name, e.g., 'S22-2' or 'NBC1-BzMe-8.0'
        self.dbrxn = dbse + '-' + str(name)
        # numerical index of reaction
        self.indx = indx
        # description line
        self.tagl = tagl
        # latex description
        self.latex = latex
        # addl comments
        self.comment = comment
        # reaction matrices, specifying reagent contributions per reaction
        self.rxnm = {}
        # qcdb.ReactionDatum objects of quantum chemical data pertaining to reaction
        self.data = {}
        # benchmark qcdb.ReactionDatum
        self.benchmark = None
        # color for plotting
        self.color = color

    def __str__(self):
        text = ''
        text += """  ==> %s Reaction <==\n\n""" % (self.name)
        text += """  Database reaction:    %s\n""" % (self.dbrxn)
        text += """  Index:                %s\n""" % (self.indx)
        text += """  LaTeX representation: %s\n""" % (self.latex)
        text += """  Tagline:              %s\n""" % (self.tagl)
        text += """  Comment:              %s\n""" % (self.comment)
        text += """  Benchmark:            %f\n""" % (self.data[self.benchmark].value)
        text += """  Color:                %s\n""" % (str(self.color))
        text += """  Reaction matrix:\n"""
        for mode, rxnm in self.rxnm.iteritems():
            text += """      %s\n""" % (mode)
            for rgt, coeff in rxnm.iteritems():
                text += """       %3d  %s\n""" % (coeff, rgt.name)
        text += """  Data:\n"""
        for label, datum in self.data.iteritems():
            text += """      %8.2f  %s\n""" % (datum.value, label)
        text += """\n"""
        return text


class Database(object):
    """

    """

    def __init__(self, dbname, pythonpath=None):
        # internal name of database
        self.dbse = None
        # OrderedDict of reactions/members
        self.hrxn = None
        # dict of reagents/geometries
        self.hrgt = None
        # dict of defined reaction subsets
        self.sset = None
        # Note that self.sset['default'] contains all the nonredundant information.
        #   Removing hrxn, hrgt etc. do not reduce the size of the object.
        #   These attributes are stored for ease of access for adding qc info, etc.

        # load database
        if pythonpath is not None:
            sys.path.insert(0, pythonpath)
        try:
            database = psiutil.import_ignorecase(dbname)
        except ImportError:
            print('\nPython module for database %s failed to load\n\n' % (dbname))
            print('\nSearch path that was tried:\n')
            print(", ".join(map(str, sys.path)))
            raise ValidationError("Python module loading problem for database " + str(dbname))

        # gross validation of database
        for item in ['dbse', 'GEOS', 'HRXN', 'ACTV', 'RXNM']:
            try:
                getattr(database, item)
            except AttributeError:
                raise ValidationError("""Database %s severely deformed with %s missing.""" % (database.__name__, item))
        for item in ['TAGL', 'BIND']:
            try:
                getattr(database, item)
            except AttributeError:
                print """Warning: Database %s possibly deformed with %s missing.\n""" % (database.__name__, item)
        self.dbse = database.dbse

        # form array of database contents to process through
        pieces = []
        for item in dir(database):
            if item in ['qcdb', 'rxn', 'dbse', 'TAGL']:
                pass
            elif item.startswith('__'):
                pass
            else:
                pieces.append(item)

        # form qcdb.Reagent objects from all defined geometries, GEOS
        oHRGT = {}
        for rgt, mol in database.GEOS.iteritems():
            mol.update_geometry()
            try:
                tagl = database.TAGL[rgt]
            except KeyError:
                tagl = None
                print """Warning: TAGL missing for reagent %s""" % (rgt)
            oHRGT[rgt] = Reagent(name=rgt, mol=mol, tagl=tagl)
        pieces.remove('GEOS')
        self.hrgt = oHRGT

        # form qcdb.Reaction objects from comprehensive reaction list, HRXN
        oHRXN = collections.OrderedDict()
        for rxn in database.HRXN:
            try:
                tagl = database.TAGL[database.dbse + '-' + str(rxn)]
            except KeyError:
                tagl = None
                print """Warning: TAGL missing for reaction %s""" % (rxn)
            try:
                elst = database.DATA['SAPT ELST ENERGY'][database.dbse + '-' + str(rxn)]
                disp = database.DATA['SAPT DISP ENERGY'][database.dbse + '-' + str(rxn)]
                color = abs(elst) / (abs(elst) + abs(disp))
            except (KeyError, AttributeError):
                color = 'black'
                print """Warning: DATA['SAPT * ENERGY'] missing for reaction %s""" % (rxn)

            oHRXN[rxn] = Reaction(name=rxn,
                                  dbse=database.dbse,
                                  indx=database.HRXN.index(rxn) + 1,
                                  color=color,
                                  tagl=tagl)
        pieces.remove('HRXN')
        self.hrxn = oHRXN

        # list and align database stoichiometry modes, ACTV* and RXNM*
        oACTV = {}
        for modactv in [item for item in pieces if item.startswith('ACTV')]:
            modrxnm = modactv.replace('ACTV', 'RXNM')
            mode = 'default' if modactv == 'ACTV' else modactv.replace('ACTV_', '')
            try:
                getattr(database, modrxnm)
            except AttributeError:
                modrxnm = 'RXNM'
            oACTV[mode] = [modactv, modrxnm]
        for item in [tmp for tmp in pieces if tmp.startswith('ACTV') or tmp.startswith('RXNM')]:
            pieces.remove(item)

        # populate reaction matrices in qcdb.Reaction objects
        for rxn in database.HRXN:
            dbrxn = database.dbse + '-' + str(rxn)
            for mode, actvrxnm in oACTV.iteritems():
                tdict = collections.OrderedDict()
                for rgt in getattr(database, actvrxnm[0])[dbrxn]:
                    tdict[oHRGT[rgt]] = getattr(database, actvrxnm[1])[dbrxn][rgt]
                oHRXN[rxn].rxnm[mode] = tdict

        # TODO neglecting case of only one BIND

        # list embedded quantum chem info per rxn, incl. BIND*
        oBIND = {}
        for arrbind in [item for item in pieces if item.startswith('BIND_')]:
            ref = arrbind.replace('BIND_', '')
            methods[ref] = Method(name=ref)
            bases[ref] = BasisSet(name=ref)
            oBIND[ref] = [methods[ref], 'default', bases[ref], (getattr(database, arrbind) is database.BIND)]
        for item in [tmp for tmp in pieces if tmp.startswith('BIND')]:
            pieces.remove(item)

        # populate data with reference values in qcdb.Reaction objects
        for rxn in database.HRXN:
            dbrxn = database.dbse + '-' + str(rxn)
            for ref, info in oBIND.iteritems():
                oHRXN[rxn].data[ref] = ReactionDatum(dbse=database.dbse,
                                                     rxn=rxn,
                                                     method=info[0],
                                                     mode=info[1],
                                                     basis=info[2],
                                                     #value=getattr(database, info[3])[dbrxn])
                                                     value=getattr(database, 'BIND_' + ref)[dbrxn])
                if info[3]:
                    oHRXN[rxn].benchmark = ref

        oSSET = {}
        fsHRXN = frozenset(database.HRXN)
        for sset in pieces:
            try:
                fssset = frozenset(getattr(database, sset))
            except TypeError:
                continue
            if fssset.issubset(fsHRXN):
                oSSET[sset] = getattr(database, sset)
        for item in oSSET.keys():
            pieces.remove(item)
        oSSET['HRXN'] = database.HRXN

        self.sset = collections.OrderedDict()
        for item in oSSET.keys():
            if item == 'HRXN_SM':
                label = 'small'
            elif item == 'HRXN_LG':
                label = 'large'
            elif item == 'HRXN_EQ':
                label = 'equilibrium'
            elif item == 'HRXN':
                label = 'default'
            elif item.startswith('HRXN_'):
                label = item.replace('HRXN_', '').lower()
            else:
                label = item.lower()

            # subsets may have different ordering from HRXN
            self.sset[label] = collections.OrderedDict()
            for rxn in oSSET[item]:
                self.sset[label][rxn] = oHRXN[rxn]

        print """Database %s: Unparsed attributes""" % (self.dbse), pieces

    def __str__(self):
        text = ''
        text += """  ==> %s Database <==\n\n""" % (self.dbse)
        text += """  Reagents:             %s\n""" % (self.hrgt.keys())
        text += """  Reactions:            %s\n""" % (self.hrxn.keys())
        text += """  Subsets:              %s\n""" % (self.sset.keys())
        text += """\n"""
        return text

    def add_ReactionDatum(self, dbse, rxn, method, mode, basis, value, units='kcal/mol', comment=None, overwrite=False):
        """Add a new quantum chemical value to *rxn* by creating a
        qcdb.ReactionDatum from same arguments as that class's
        object-less constructor. *rxn* may be actual Reaction.name
        or Reaction.indx.

        """
        if (self.dbse == dbse):
            if rxn in self.hrxn.keys():
                rxnname = rxn  # rxn is proper reaction name
            else:
                try:
                    if (rxn + 1 > 0) and (rxn == self.hrxn.items()[rxn - 1][1].indx):
                        rxnname = self.hrxn.items()[rxn - 1][1].name  # rxn is reaction index (maybe dangerous?)
                except (TypeError, IndexError):
                    raise ValidationError("""Inconsistent to add ReactionDatum for %s to database %s with reactions %s.""" %
                        (dbse + '-' + str(rxn), self.dbse, self.hrxn.keys()))
            label = '-'.join([method, mode, basis])
            if overwrite or (label not in self.hrxn[rxnname].data.keys()):
                self.hrxn[rxnname].data[label] = ReactionDatum.library_modelchem(dbse=dbse, rxn=rxnname,
                    method=method, mode=mode, basis=basis,
                    value=value, units=units, comment=comment)
            else:
                raise ValidationError("""ReactionDatum %s already present in Database.""" % (label))
        else:
            raise ValidationError("""Inconsistent to add ReactionDatum for %s to database %s.""" %
                (dbse + '-' + str(rxn), self.dbse))

    def add_Subset(self, name, func):
        """Define a new subset labeled *name* by providing a function
        *func* that filters *self.hrxn*.

        """
        label = name.lower()
        try:
            lsslist = [rxn for rxn in self.sset['default'].keys() if rxn in func(self)]
        except TypeError, e:
            raise ValidationError("""Function %s did not return list: %s.""" % (func.__name__, str(e)))

        self.sset[label] = collections.OrderedDict()
        for rxn in lsslist:
            self.sset[label][rxn] = self.hrxn[rxn]
        print """Database %s: Subset %s formed: %s""" % (self.dbse, label, self.sset[label].keys())

    def compute_errors(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False):
        """For full database or subset *sset*, computes raw reaction
        errors between *modelchem* and *benchmark* model chemistries.
        Returns error if model chemistries are missing for any reaction in
        subset unless *failoninc* set to False, whereupon returns partial.
        Returns dictionary of reaction labels and error forms.

        """
        if isinstance(sset, basestring):
            # sset is normal subset name 'MX' corresponding to HRXN_MX or MX array in database module
            try:
                lsset = self.sset[sset.lower()]
            except KeyError, e:
                #raise ValidationError("""Subset named %s not available""" % (str(e)))
                lsset = collections.OrderedDict()
        else:
            if callable(sset):
                # sset is function that will generate subset of HRXN from sset(self)
                lsslist = [rxn for rxn in self.sset['default'].keys() if rxn in sset(self)]
            else:
                # sset is array containing reactions
                lsslist = [rxn for rxn in self.sset['default'].keys() if rxn in sset]
            # assemble dict of qcdb.Reaction objects from array of reaction names
            lsset = collections.OrderedDict()
            for rxn in lsslist:
                lsset[rxn] = self.hrxn[rxn]

        err = {}
        for rxn, oRxn in lsset.iteritems():
            lbench = oRxn.benchmark if benchmark == 'default' else benchmark
            try:
                mcLesser = oRxn.data[modelchem].value
                mcGreater = oRxn.data[lbench].value
            except KeyError, e:
                if failoninc:
                    raise ValidationError("""Reaction %s missing datum %s.""" % (str(rxn), str(e)))
                else:
                    continue

            err[rxn] = [mcLesser - mcGreater,
                        (mcLesser - mcGreater) / abs(mcGreater),
                        (mcLesser - mcGreater) / abs(mcGreater)]  # TODO define BER
            if verbose:
                print """p = %6.2f, pe = %6.1f%%, bpe = %6.1f%% reaction %s.""" % \
                    (err[rxn][0], 100 * err[rxn][1], 100 * err[rxn][2], str(rxn))
        return err

    def compute_statistics(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False, returnindiv=False):
        """For full database or subset *sset*, computes many error
        statistics between single *modelchem* and *benchmark* model
        chemistries. Returns error if model chemistries are missing
        for any reaction in subset unless *failoninc* set to False,
        whereupon returns partial statistics. Returns dictionary of
        statistics labels and values.

        """
        err = self.compute_errors(modelchem, benchmark=benchmark, sset=sset, failoninc=failoninc, verbose=verbose)
        if len(err) == 0:
            error = initialize_errors()
            if verbose:
                print """Warning: nothing to compute."""
        else:
            Nrxn = float(len(err))
            error = collections.OrderedDict()
            # linear (absolute) error
            linear = [val[0] for val in err.values()]
            error['maxe'] = max(linear, key=lambda x: abs(x))
            error['mine'] = min(linear, key=lambda x: abs(x))
            error['me'] = sum(linear) / Nrxn
            error['mae'] = sum(map(abs, linear)) / Nrxn
            error['rmse'] = math.sqrt(sum(map(lambda x: x ** 2, linear)) / Nrxn)
            error['stde'] = math.sqrt((sum(map(lambda x: x ** 2, linear)) - (sum(linear) ** 2) / Nrxn) / Nrxn)
            # fractional (relative) error
            relative = [val[1] for val in err.values()]
            error['maxpe'] = max(relative, key=lambda x: abs(x))
            error['minpe'] = min(relative, key=lambda x: abs(x))
            error['mpe'] = sum(relative) / Nrxn
            error['mape'] = sum(map(abs, relative)) / Nrxn
            error['rmspe'] = math.sqrt(sum(map(lambda x: x ** 2, relative)) / Nrxn)
            error['stdpe'] = math.sqrt((sum(map(lambda x: x ** 2, relative)) - (sum(relative) ** 2) / Nrxn) / Nrxn)
            # balanced (relative) error
            balanced = [val[2] for val in err.values()]
            error['maxpbe'] = max(balanced, key=lambda x: abs(x))
            error['minpbe'] = min(balanced, key=lambda x: abs(x))
            error['mpbe'] = sum(balanced) / Nrxn
            error['mapbe'] = sum(map(abs, balanced)) / Nrxn
            error['rmspbe'] = math.sqrt(sum(map(lambda x: x ** 2, balanced)) / Nrxn)
            error['stdpbe'] = math.sqrt((sum(map(lambda x: x ** 2, balanced)) - (sum(balanced) ** 2) / Nrxn) / Nrxn)
            if verbose:
                print """%d systems in %s for %s vs. %s, subset %s.\n%s""" % \
                    (len(err), self.dbse, modelchem, benchmark, sset, format_errors(error, mode=2))
        if returnindiv:
            return error, err
        else:
            return error

    def load_qcdata(self, modname, funcname, pythonpath=None, failoninc=True):
        """Loads qcdb.ReactionDatums from module *modname* function
        *funcname*. Module search path can be prepended with *pythonpath*.

        """
        if pythonpath is not None:
            sys.path.insert(0, pythonpath)
        else:
            sys.path.append(os.path.dirname(__file__) + '/../data')
        try:
            datamodule = __import__(modname)
        except ImportError:
            if not failoninc:
                print """%s data unavailable for database %s.\n""" % (modname, self.dbse)
                return
            else:
                print '\nPython module for database data %s failed to load\n\n' % (modname)
                print '\nSearch path that was tried:\n'
                print ", ".join(map(str, sys.path))
                raise ValidationError("Python module loading problem for database data " + str(modname))
        try:
            getattr(datamodule, funcname)(self)
        except AttributeError:
            if not failoninc:
                print """%s %s data unavailable for database %s.\n""" % (modname, funcname, self.dbse)
                return
            else:
                raise ValidationError("Python module missing function %s for loading data " % (str(funcname)))

        print """Database %s: %s %s results loaded""" % (self.dbse, modname, funcname)

    def load_qcdata_byproject(self, project, pythonpath=None):
        """Loads qcdb.ReactionDatums from standard location for *project*
        :module dbse_project and function load_project. Module search path
        can be prepended with *pythonpath*.

        """
        mod = self.dbse + '_' + project
        func = 'load_' + project
        self.load_qcdata(modname=mod, funcname=func, pythonpath=pythonpath)

    def load_subsets(self, modname='subsetgenerator', pythonpath=None):
        """Loads subsets from all functions in module *modname*.

        """
        if pythonpath is not None:
            sys.path.insert(0, pythonpath)
        else:
            sys.path.append(os.path.dirname(__file__))
        try:
            ssmod = __import__(modname)
        except ImportError:
            print '\nPython module for database data %s failed to load\n\n' % (modname)
            print '\nSearch path that was tried:\n'
            print ", ".join(map(str, sys.path))
            raise ValidationError("Python module loading problem for database subset generator " + str(modname))

        for func in dir(ssmod):
            if callable(getattr(ssmod, func)):
                self.add_Subset(getattr(ssmod, func).__doc__, getattr(ssmod, func))

        print """Database %s: Defined subsets loaded""" % (self.dbse)

    #def analyze_modelchem(self, modelchem, benchmark='default', failoninc=True, verbose=False):
    #    """Compute and print error statistics for *modelchem* versus
    #    *benchmark* for all available subsets and return dictionary of same.

    #    """
    #    errors = collections.OrderedDict()
    #    for ss in self.sset.keys():
    #        errors[ss] = self.compute_statistics(modelchem, benchmark=benchmark, sset=ss, failoninc=failoninc, verbose=verbose)
    #    print """\n  ==> %s %s Errors <==""" % (self.dbse, modelchem)
    #    print """%20s        %5s  %4s   %6s %6s    %6s""" % \
    #        ('', 'ME', 'STDE', 'MAE', 'MA%E', 'MA%BE')
    #    for ss in errors.keys():
    #        if any(errors[ss].values()):
    #            print """%20s    %42s""" % (ss, format_errors(errors[ss]))
    #    return errors

    def analyze_modelchems(self, modelchem, benchmark='default', failoninc=True, verbose=False):
        """Compute and print error statistics for each model chemistry in
        array *modelchem* versus *benchmark* for all available subsets and
        return dictionary of same.

        """
        pre, suf, mid = string_contrast(modelchem)
        errors = collections.OrderedDict()
        for ss in self.sset.keys():
            errors[ss] = collections.OrderedDict()
            for mc in modelchem:
                errors[ss][mc] = self.compute_statistics(mc, benchmark=benchmark, sset=ss, failoninc=failoninc, verbose=verbose)
        print """\n  ==> %s %s[]%s Errors <==""" % (self.dbse, pre, suf)
        print """%20s        %5s  %4s   %6s %6s    %6s""" % \
            ('', 'ME', 'STDE', 'MAE', 'MA%E', 'MA%BE')
        for ss in self.sset.keys():
            if any([any(errors[ss][mc].values()) for mc in modelchem]):
                print """   => %s <= """ % (ss)
                for mc in modelchem:
                    print """%20s    %42s""" % (mid[modelchem.index(mc)], format_errors(errors[ss][mc]))
        return errors

    def plot_modelchems(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False, color='sapt', xlimit=4.0):
        """Computes individual errors and summary statistics for each model
        chemistry in array *modelchem* versus *benchmark* over subset *sset*.
        Thread *color* can be 'rgb' for old coloring, a color name or 'sapt'
        for spectrum coloring. Prepares thread diagram instructions and
        either executes them if matplotlib available (Canopy) or prints them.

        """
        pre, suf, mid = string_contrast(modelchem)
        # compute errors
        errors = collections.OrderedDict()
        indiv = collections.OrderedDict()
        for mc in modelchem:
            errors[mc], indiv[mc] = self.compute_statistics(mc, benchmark=benchmark,
                sset=sset, failoninc=failoninc, verbose=verbose, returnindiv=True)
        #print 'ERRORS', errors
        # repackage
        dbdat = []
        for rxn in self.sset[sset].keys():
            data = []
            for mc in modelchem:
                try:
                    data.append(indiv[mc][rxn][0])
                except KeyError, e:
                    if failoninc:
                        raise e
                    else:
                        data.append(None)
            dbdat.append({'sys': str(rxn), 'color': self.hrxn[rxn].color, 'data': data})
        #dbdat = [{'sys': str(rxn), 'color': self.hrxn[rxn].color,
        #    'data': [indiv[mc][rxn][0] for mc in modelchem]} for rxn in self.sset[sset].keys()]
        title = self.dbse + ' ' + pre + '[]' + suf
        mae = [errors[mc]['mae'] for mc in modelchem]
        mapbe = [100 * errors[mc]['mapbe'] for mc in modelchem]
        # generate matplotlib instructions and call or print
        try:
            import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print """mpl.thread(%s,\n    color='%s',\n    title='%s',\n    labels=%s,\n    mae=%s,\n    mape=%s\n    xlimit=%s)\n\n""" % \
                (dbdat, color, title, mid, mae, mapbe, str(xlimit))
        else:
            # if running from Canopy, call mpl directly
            mpl.thread(dbdat, color=color, title=title, labels=mid, mae=mae, mape=mapbe, xlimit=xlimit)
            print """mpl.thread(%s,\n    color='%s',\n    title='%s',\n    labels=%s,\n    mae=%s,\n    mape=%s\n    xlimit=%s)\n\n""" % \
                (dbdat, color, title, mid, mae, mapbe, str(xlimit))

    def plot_bars(self, modelchem, benchmark='default', sset=['default', 'hb', 'mx', 'dd'], failoninc=True, verbose=False):
        """Prepares 'grey bars' diagram for each model chemistry in array
        *modelchem* versus *benchmark* over all four databases. A wide bar
        is plotted with three smaller bars, corresponding to the 'mae'
        summary statistic of the four subsets in *sset*. Prepares bars
        diagram instructions and either executes them if matplotlib
        available (Canopy) or prints them.

        """
        # compute errors
        errors = {}
        for mc in modelchem:
            if mc is not None:
                errors[mc] = {}
                for ss in sset:
                    errors[mc][ss] = self.compute_statistics(mc, benchmark=benchmark, sset=ss,
                        failoninc=failoninc, verbose=verbose, returnindiv=False)
        # repackage
        pre, suf, mid = string_contrast(modelchem)
        #dbdat = [{'mc': mid[modelchem.index(mc)], 'data': [errors[mc][ss]['DB4']['mae'] for ss in sset]} for mc in modelchem]
        dbdat = []
        for mc in modelchem:
            if mc is None:
                dbdat.append(None)
            else:
                dbdat.append({'mc': mid[modelchem.index(mc)], 'data': [errors[mc][ss]['mae'] for ss in sset]})
        title = self.dbse + ' ' + pre + '[]' + suf
        # generate matplotlib instructions and call or print
        try:
            import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print """mpl.bar(%s,\n    title='%s')\n\n""" % (dbdat, title)
        else:
            # if running from Canopy, call mpl directly
            mpl.bar(dbdat, title=title)


class FourDatabases(object):
    """

    """
    def __init__(self):
        sys.path.append(os.path.dirname(__file__) + '/../databases')
        # S22 database
        self.s22 = Database('S22')
        # NBC10 database
        self.nbc1 = Database('NBC10')
        # HBC6 database
        self.hbc1 = Database('HBC6')
        # HSG database
        self.hsg = Database('HSG')
        # subset assembly pattern
        self.sset = collections.OrderedDict()
        # assembly pattern for transspecies modelchems
        self.mc = {}

        # load up data and definitions
        #self.load_qcdata_byproject('dft')
        self.load_qcdata_byproject('pt2')
        #self.load_qcdata_byproject('dhdft')
        self.load_subsets()
        self.define_supersubsets()
        self.define_supermodelchems()

    def load_qcdata_byproject(self, project, pythonpath=None):
        """Load data for *project* from standard location, modifiable by
        *pythonpath*.

        """
        self.s22.load_qcdata_byproject(project, pythonpath=pythonpath)
        self.nbc1.load_qcdata_byproject(project, pythonpath=pythonpath)
        self.hbc1.load_qcdata_byproject(project, pythonpath=pythonpath)
        self.hsg.load_qcdata_byproject(project, pythonpath=pythonpath)

    def load_subsets(self):
        """Load subsets from standard generators.

        """
        self.s22.load_subsets()
        self.nbc1.load_subsets()
        self.hbc1.load_subsets()
        self.hsg.load_subsets()

    def define_supersubsets(self):
        """

        """
        self.sset['default'] = ['default', 'default', 'default', 'default']
        self.sset['tt'] = ['default', 'default', 'default', 'default']
        self.sset['hb'] = ['hb', None, 'default', 'hb']
        self.sset['mx'] = ['mx', 'mx', None, 'mx']
        self.sset['dd'] = ['dd', 'dd', None, 'dd']
        self.sset['mxdd'] = ['mxdd', 'default', None, 'mxdd']
        self.sset['pp'] = ['mxddpp', 'mxddpp', None, None]
        self.sset['np'] = ['mxddnp', 'mxddnp', None, 'mxdd']
        self.sset['tt-5min'] = ['default', '5min', '5min', 'default']
        self.sset['hb-5min'] = ['hb', None, '5min', 'hb']
        self.sset['mx-5min'] = ['mx', 'mx-5min', None, 'mx']
        self.sset['dd-5min'] = ['dd', 'dd-5min', None, 'dd']
        self.sset['mxdd-5min'] = ['mxdd', '5min', None, 'mxdd']
        self.sset['pp-5min'] = ['mxddpp', 'mxddpp-5min', None, None]
        self.sset['np-5min'] = ['mxddnp', 'mxddnp-5min', None, 'mxdd']

    def define_supermodelchems(self):
        """

        """
        self.mc['C2011BENCH'] = ['S22A', 'NBC100', 'HBC60', 'HSG0']

        self.mc['CCSD-CP-adz'] = ['CCSD-CP-adz', 'CCSD-CP-hadz', 'CCSD-CP-adz', 'CCSD-CP-hadz']
        self.mc['CCSD-CP-atz'] = ['CCSD-CP-atz', 'CCSD-CP-hatz', 'CCSD-CP-atz', 'CCSD-CP-hatz']
        self.mc['CCSD-CP-adtz'] = ['CCSD-CP-adtz', 'CCSD-CP-hadtz', 'CCSD-CP-adtz', 'CCSD-CP-hadtz']
        self.mc['CCSD-CP-adtzadz'] = ['CCSD-CP-adtzadz', 'CCSD-CP-adtzhadz', 'CCSD-CP-adtzadz', 'CCSD-CP-adtzhadz']
        self.mc['CCSD-CP-atzadz'] = ['CCSD-CP-atzadz', 'CCSD-CP-atzhadz', 'CCSD-CP-atzadz', 'CCSD-CP-atzhadz']
        self.mc['CCSD-CP-atqzadz'] = ['CCSD-CP-atqzadz', 'CCSD-CP-atqzhadz', 'CCSD-CP-atqzadz', 'CCSD-CP-atqzhadz']
        self.mc['CCSD-CP-atzadtz'] = ['CCSD-CP-atzadtz', 'CCSD-CP-atzhadtz', 'CCSD-CP-atzadtz', 'CCSD-CP-atzhadtz']
        self.mc['CCSD-CP-atqzadtz'] = ['CCSD-CP-atqzadtz', 'CCSD-CP-atqzhadtz', 'CCSD-CP-atqzadtz', 'CCSD-CP-atqzhadtz']
        self.mc['CCSD-CP-atqzatz'] = ['CCSD-CP-atqzatz', 'CCSD-CP-atqzhatz', 'CCSD-CP-atqzatz', 'CCSD-CP-atqzhatz']

        self.mc['SCSCCSD-CP-adz'] = ['SCSCCSD-CP-adz', 'SCSCCSD-CP-hadz', 'SCSCCSD-CP-adz', 'SCSCCSD-CP-hadz']
        self.mc['SCSCCSD-CP-atz'] = ['SCSCCSD-CP-atz', 'SCSCCSD-CP-hatz', 'SCSCCSD-CP-atz', 'SCSCCSD-CP-hatz']
        self.mc['SCSCCSD-CP-adtz'] = ['SCSCCSD-CP-adtz', 'SCSCCSD-CP-hadtz', 'SCSCCSD-CP-adtz', 'SCSCCSD-CP-hadtz']
        self.mc['SCSCCSD-CP-adtzadz'] = ['SCSCCSD-CP-adtzadz', 'SCSCCSD-CP-adtzhadz', 'SCSCCSD-CP-adtzadz', 'SCSCCSD-CP-adtzhadz']
        self.mc['SCSCCSD-CP-atzadz'] = ['SCSCCSD-CP-atzadz', 'SCSCCSD-CP-atzhadz', 'SCSCCSD-CP-atzadz', 'SCSCCSD-CP-atzhadz']
        self.mc['SCSCCSD-CP-atqzadz'] = ['SCSCCSD-CP-atqzadz', 'SCSCCSD-CP-atqzhadz', 'SCSCCSD-CP-atqzadz', 'SCSCCSD-CP-atqzhadz']
        self.mc['SCSCCSD-CP-atzadtz'] = ['SCSCCSD-CP-atzadtz', 'SCSCCSD-CP-atzhadtz', 'SCSCCSD-CP-atzadtz', 'SCSCCSD-CP-atzhadtz']
        self.mc['SCSCCSD-CP-atqzadtz'] = ['SCSCCSD-CP-atqzadtz', 'SCSCCSD-CP-atqzhadtz', 'SCSCCSD-CP-atqzadtz', 'SCSCCSD-CP-atqzhadtz']
        self.mc['SCSCCSD-CP-atqzatz'] = ['SCSCCSD-CP-atqzatz', 'SCSCCSD-CP-atqzhatz', 'SCSCCSD-CP-atqzatz', 'SCSCCSD-CP-atqzhatz']

        self.mc['SCSMICCSD-CP-adz'] = ['SCSMICCSD-CP-adz', 'SCSMICCSD-CP-hadz', 'SCSMICCSD-CP-adz', 'SCSMICCSD-CP-hadz']
        self.mc['SCSMICCSD-CP-atz'] = ['SCSMICCSD-CP-atz', 'SCSMICCSD-CP-hatz', 'SCSMICCSD-CP-atz', 'SCSMICCSD-CP-hatz']
        self.mc['SCSMICCSD-CP-adtz'] = ['SCSMICCSD-CP-adtz', 'SCSMICCSD-CP-hadtz', 'SCSMICCSD-CP-adtz', 'SCSMICCSD-CP-hadtz']
        self.mc['SCSMICCSD-CP-adtzadz'] = ['SCSMICCSD-CP-adtzadz', 'SCSMICCSD-CP-adtzhadz', 'SCSMICCSD-CP-adtzadz', 'SCSMICCSD-CP-adtzhadz']
        self.mc['SCSMICCSD-CP-atzadz'] = ['SCSMICCSD-CP-atzadz', 'SCSMICCSD-CP-atzhadz', 'SCSMICCSD-CP-atzadz', 'SCSMICCSD-CP-atzhadz']
        self.mc['SCSMICCSD-CP-atqzadz'] = ['SCSMICCSD-CP-atqzadz', 'SCSMICCSD-CP-atqzhadz', 'SCSMICCSD-CP-atqzadz', 'SCSMICCSD-CP-atqzhadz']
        self.mc['SCSMICCSD-CP-atzadtz'] = ['SCSMICCSD-CP-atzadtz', 'SCSMICCSD-CP-atzhadtz', 'SCSMICCSD-CP-atzadtz', 'SCSMICCSD-CP-atzhadtz']
        self.mc['SCSMICCSD-CP-atqzadtz'] = ['SCSMICCSD-CP-atqzadtz', 'SCSMICCSD-CP-atqzhadtz', 'SCSMICCSD-CP-atqzadtz', 'SCSMICCSD-CP-atqzhadtz']
        self.mc['SCSMICCSD-CP-atqzatz'] = ['SCSMICCSD-CP-atqzatz', 'SCSMICCSD-CP-atqzhatz', 'SCSMICCSD-CP-atqzatz', 'SCSMICCSD-CP-atqzhatz']

        self.mc['CCSDT-CP-adz'] = ['CCSDT-CP-adz', 'CCSDT-CP-hadz', 'CCSDT-CP-adz', 'CCSDT-CP-hadz']
        self.mc['CCSDT-CP-atz'] = ['CCSDT-CP-atz', 'CCSDT-CP-hatz', 'CCSDT-CP-atz', 'CCSDT-CP-hatz']
        self.mc['CCSDT-CP-adtz'] = ['CCSDT-CP-adtz', 'CCSDT-CP-hadtz', 'CCSDT-CP-adtz', 'CCSDT-CP-hadtz']
        self.mc['CCSDT-CP-adtzadz'] = ['CCSDT-CP-adtzadz', 'CCSDT-CP-adtzhadz', 'CCSDT-CP-adtzadz', 'CCSDT-CP-adtzhadz']
        self.mc['CCSDT-CP-atzadz'] = ['CCSDT-CP-atzadz', 'CCSDT-CP-atzhadz', 'CCSDT-CP-atzadz', 'CCSDT-CP-atzhadz']
        self.mc['CCSDT-CP-atqzadz'] = ['CCSDT-CP-atqzadz', 'CCSDT-CP-atqzhadz', 'CCSDT-CP-atqzadz', 'CCSDT-CP-atqzhadz']
        self.mc['CCSDT-CP-atzadtz'] = ['CCSDT-CP-atzadtz', 'CCSDT-CP-atzhadtz', 'CCSDT-CP-atzadtz', 'CCSDT-CP-atzhadtz']
        self.mc['CCSDT-CP-atqzadtz'] = ['CCSDT-CP-atqzadtz', 'CCSDT-CP-atqzhadtz', 'CCSDT-CP-atqzadtz', 'CCSDT-CP-atqzhadtz']
        self.mc['CCSDT-CP-atqzatz'] = ['CCSDT-CP-atqzatz', 'CCSDT-CP-atqzhatz', 'CCSDT-CP-atqzatz', 'CCSDT-CP-atqzhatz']

    def analyze_modelchems(self, modelchem, benchmark='default', failoninc=True, verbose=False):
        """Print out nicely formatted summary statistics for every model
        chemistry in array *modelchem* versus *benchmark* for every
        registered subset.

        """
        # compute errors
        errors = {}
        for mc in modelchem:
            errors[mc] = {}
            for ss in self.sset.keys():
                errors[mc][ss] = self.compute_statistics(mc, benchmark=benchmark, sset=ss,
                    failoninc=failoninc, verbose=verbose, returnindiv=False)
        # present errors
        pre, suf, mid = string_contrast(modelchem)
        print """\n  ==> %s %s[]%s Errors <==""" % ('DB4', pre, suf)
        print """%20s    %42s%42s%42s%42s%42s""" % \
            ('', '=> DB4 <=', '=> S22 <=', '=> NBC1 <=', '=> HBC1 <=', '=> HSG <=')
        print """%20s        %5s  %4s   %6s %6s    %6s""" % \
            ('', 'ME', 'STDE', 'MAE', 'MA%E', 'MA%BE')
        for ss in self.sset.keys():
            print """   => %s <= """ % (ss)
            for mc in modelchem:
                tmpdb = errors[mc][ss]
                print """%20s    %42s%42s%42s%42s%42s""" % (mid[modelchem.index(mc)],
                    format_errors(errors[mc][ss]['DB4']),
                    '' if tmpdb['S22'] is None else format_errors(tmpdb['S22']),
                    '' if tmpdb['NBC1'] is None else format_errors(tmpdb['NBC1']),
                    '' if tmpdb['HBC1'] is None else format_errors(tmpdb['HBC1']),
                    '' if tmpdb['HSG'] is None else format_errors(tmpdb['HSG']))

    def plot_bars(self, modelchem, benchmark='default', sset=['tt-5min', 'hb-5min', 'mx-5min', 'dd-5min'], failoninc=True, verbose=False):
        """Prepares 'grey bars' diagram for each model chemistry in array
        *modelchem* versus *benchmark* over all four databases. A wide bar
        is plotted with three smaller bars, corresponding to the 'mae'
        summary statistic of the four subsets in *sset*. Prepares bars
        diagram instructions and either executes them if matplotlib
        available (Canopy) or prints them.

        """
        # compute errors
        errors = {}
        for mc in modelchem:
            if mc is not None:
                errors[mc] = {}
                for ss in sset:
                    errors[mc][ss] = self.compute_statistics(mc, benchmark=benchmark, sset=ss,
                        failoninc=failoninc, verbose=verbose, returnindiv=False)
        # repackage
        pre, suf, mid = string_contrast(modelchem)
        #dbdat = [{'mc': mid[modelchem.index(mc)], 'data': [errors[mc][ss]['DB4']['mae'] for ss in sset]} for mc in modelchem]
        dbdat = []
        for mc in modelchem:
            if mc is None:
                dbdat.append(None)
            else:
                dbdat.append({'mc': mid[modelchem.index(mc)], 'data': [errors[mc][ss]['DB4']['mae'] for ss in sset]})
        title = '4DB ' + pre + '[]' + suf
        # generate matplotlib instructions and call or print
        try:
            import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print """mpl.bar(%s,\n    title='%s')\n\n""" % (dbdat, title)
        else:
            # if running from Canopy, call mpl directly
            mpl.bar(dbdat, title=title)

    def compute_statistics(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False, returnindiv=False):
        """Computes summary statistics and, if *returnindiv* True,
        individual errors for single model chemistry *modelchem* versus
        *benchmark* over subset *sset* over all four databases.
        Particularly, imposes cross-database definitions for sset and
        modelchem.

        """
        pre, suf, mid = string_contrast(modelchem)
        if modelchem in self.mc.keys():
            lmc = self.mc[modelchem]
        else:
            lmc = [modelchem, modelchem, modelchem, modelchem]
        if benchmark in self.mc.keys():
            lbm = self.mc[benchmark]
        else:
            lbm = [benchmark, benchmark, benchmark, benchmark]
        if sset in self.sset.keys():
            lss = self.sset[sset]
        else:
            lss = [sset, sset, sset, sset]

        errors = collections.OrderedDict()
        indiv = collections.OrderedDict()
        errors['S22'], indiv['S22'] = (None, None) if lss[0] is None else self.s22.compute_statistics(lmc[0], sset=lss[0],
            benchmark=lbm[0], failoninc=failoninc, verbose=verbose, returnindiv=True)
        errors['NBC1'], indiv['NBC1'] = (None, None) if lss[1] is None else self.nbc1.compute_statistics(lmc[1], sset=lss[1],
            benchmark=lbm[1], failoninc=failoninc, verbose=verbose, returnindiv=True)
        errors['HBC1'], indiv['HBC1'] = (None, None) if lss[2] is None else self.hbc1.compute_statistics(lmc[2], sset=lss[2],
            benchmark=lbm[2], failoninc=failoninc, verbose=verbose, returnindiv=True)
        errors['HSG'], indiv['HSG'] = (None, None) if lss[3] is None else self.hsg.compute_statistics(lmc[3], sset=lss[3],
            benchmark=lbm[3], failoninc=failoninc, verbose=verbose, returnindiv=True)

        args = []
        if lss[0] is not None:
            args.append(errors['S22'])
        if lss[1] is not None:
            args.append(errors['NBC1'])
        if lss[2] is not None:
            args.append(errors['HBC1'])
        if lss[3] is not None:
            args.append(errors['HSG'])
        errors['DB4'] = average_errors(*args)

        if returnindiv:
            return errors, indiv
        else:
            return errors

    def plot_modelchems(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False, color='sapt', xlimit=4.0):
        """Computes individual errors and summary statistics for each
        model chemistry in array *modelchem* versus *benchmark* over
        subset *sset* over all four databases. Thread *color* can be 'rgb'
        for old coloring, a color name or 'sapt' for spectrum coloring.
        Prepares thread diagram instructions and either executes them if
        matplotlib available (Canopy) or prints them.

        """
        # compute errors
        errors = {}
        indiv = {}
        for mc in modelchem:
            errors[mc], indiv[mc] = self.compute_statistics(mc, benchmark=benchmark, sset=sset,
                failoninc=failoninc, verbose=verbose, returnindiv=True)
        # repackage
        dbdat = []
        for db in indiv[modelchem[0]].keys():
            if indiv[modelchem[0]][db] is not None:
                for rxn, orxn in indiv[modelchem[0]][db].items():
                    if db == 'S22':
                        lcolor = self.s22.hrxn[rxn].color
                    elif db == 'NBC1':
                        lcolor = self.nbc1.hrxn[rxn].color
                    elif db == 'HBC1':
                        lcolor = self.hbc1.hrxn[rxn].color
                    elif db == 'HSG':
                        lcolor = self.hsg.hrxn[rxn].color
                    dbdat.append({'sys': str(rxn), 'color': lcolor,
                        'data': [indiv[mc][db][rxn][0] for mc in modelchem]})
        pre, suf, mid = string_contrast(modelchem)
        title = '4DB-' + sset + ' ' + pre + '[]' + suf
        mae = [errors[mc]['DB4']['mae'] for mc in modelchem]
        mapbe = [100 * errors[mc]['DB4']['mapbe'] for mc in modelchem]
        # generate matplotlib instructions and call or print
        try:
            import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print """mpl.thread(%s,\n    color='%s',\n    title='%s',\n    labels=%s,\n    mae=%s,\n    mape=%s\n    xlimit=%s)\n\n""" % \
                (dbdat, color, title, mid, mae, mapbe, str(xlimit))
        else:
            # if running from Canopy, call mpl directly
            mpl.thread(dbdat, color=color, title=title, labels=mid, mae=mae, mape=mapbe, xlimit=xlimit)

    def plot_flat(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False, color='sapt', xlimit=4.0, view=True):
        """Computes individual errors and summary statistics for single
        model chemistry *modelchem* versus *benchmark* over
        subset *sset* over all four databases. Thread *color* can be 'rgb'
        for old coloring, a color name or 'sapt' for spectrum coloring.
        Prepares flat diagram instructions and either executes them if
        matplotlib available (Canopy) or prints them.

        """
        # compute errors
        errors = {}
        indiv = {}
        mc = modelchem
        errors[mc], indiv[mc] = self.compute_statistics(mc, benchmark=benchmark, sset=sset,
            failoninc=failoninc, verbose=verbose, returnindiv=True)
        # repackage
        dbdat = []
        for db in indiv[mc].keys():
            if indiv[mc][db] is not None:
                for rxn, orxn in indiv[mc][db].items():
                    if db == 'S22':
                        lcolor = self.s22.hrxn[rxn].color
                    elif db == 'NBC1':
                        lcolor = self.nbc1.hrxn[rxn].color
                    elif db == 'HBC1':
                        lcolor = self.hbc1.hrxn[rxn].color
                    elif db == 'HSG':
                        lcolor = self.hsg.hrxn[rxn].color
                    dbdat.append({'sys': str(rxn), 'color': lcolor,
                        'data': [indiv[mc][db][rxn][0]]})
        pre, suf, mid = string_contrast(mc)
        title = '4DB-' + sset + ' ' + pre + '[]' + suf
        mae = errors[mc]['DB4']['mae']
        mapbe = 100 * errors[mc]['DB4']['mapbe']
        mapbe = None
        # generate matplotlib instructions and call or print
        try:
            import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print """mpl.flat(%s,\n    color='%s',\n    title='%s',\n    mae=%s,\n    mape=%s,\n    xlimit=%s,\n    view=%s)\n\n""" % \
                (dbdat, color, mc, mae, mapbe, xlimit, view)
        else:
            # if running from Canopy, call mpl directly
            mpl.flat(dbdat, color=color, title=mc, mae=mae, mape=mapbe, xlimit=xlimit, view=view)

    def plot_all_flats(self):
        """Generate pieces for inclusion into tables for PT2 paper."""
        for mc in sorted(self.s22.hrxn[2].data.keys()):
            if mc not in ['S220', 'S22A', 'S22B']:
                self.plot_flat(mc, sset='tt-5min', xlimit=4.0, view=False)

    def make_pt2_Figure_3(self):
        """Plot all the graphics needed for the calendar grey bars plot
        in Fig. 3 of PT2.

        """
        # Fig. bars (a)
        self.plot_bars(['MP2-CP-dz', 'MP2-CP-jadz', 'MP2-CP-hadz', 'MP2-CP-adz',
            'MP2-CP-tz', 'MP2-CP-matz', 'MP2-CP-jatz', 'MP2-CP-hatz', 'MP2-CP-atz',
            'MP2-CP-dtz', 'MP2-CP-jadtz', 'MP2-CP-hadtz', 'MP2-CP-adtz',
            'MP2-CP-qz', 'MP2-CP-aaqz', 'MP2-CP-maqz', 'MP2-CP-jaqz', 'MP2-CP-haqz', 'MP2-CP-aqz',
            'MP2-CP-tqz', 'MP2-CP-matqz', 'MP2-CP-jatqz', 'MP2-CP-hatqz', 'MP2-CP-atqz',
            'MP2-CP-a5z', 'MP2-CP-aq5z'])
        self.plot_bars(['SCSMP2-CP-dz', 'SCSMP2-CP-jadz', 'SCSMP2-CP-hadz', 'SCSMP2-CP-adz',
            'SCSMP2-CP-tz', 'SCSMP2-CP-matz', 'SCSMP2-CP-jatz', 'SCSMP2-CP-hatz', 'SCSMP2-CP-atz',
            'SCSMP2-CP-dtz', 'SCSMP2-CP-jadtz', 'SCSMP2-CP-hadtz', 'SCSMP2-CP-adtz',
            'SCSMP2-CP-qz', 'SCSMP2-CP-aaqz', 'SCSMP2-CP-maqz', 'SCSMP2-CP-jaqz', 'SCSMP2-CP-haqz', 'SCSMP2-CP-aqz',
            'SCSMP2-CP-tqz', 'SCSMP2-CP-matqz', 'SCSMP2-CP-jatqz', 'SCSMP2-CP-hatqz', 'SCSMP2-CP-atqz',
            'SCSMP2-CP-a5z', 'SCSMP2-CP-aq5z'])
        self.plot_bars(['SCSNMP2-CP-dz', 'SCSNMP2-CP-jadz', 'SCSNMP2-CP-hadz', 'SCSNMP2-CP-adz',
            'SCSNMP2-CP-tz', 'SCSNMP2-CP-matz', 'SCSNMP2-CP-jatz', 'SCSNMP2-CP-hatz', 'SCSNMP2-CP-atz',
            'SCSNMP2-CP-dtz', 'SCSNMP2-CP-jadtz', 'SCSNMP2-CP-hadtz', 'SCSNMP2-CP-adtz',
            'SCSNMP2-CP-qz', 'SCSNMP2-CP-aaqz', 'SCSNMP2-CP-maqz', 'SCSNMP2-CP-jaqz', 'SCSNMP2-CP-haqz', 'SCSNMP2-CP-aqz',
            'SCSNMP2-CP-tqz', 'SCSNMP2-CP-matqz', 'SCSNMP2-CP-jatqz', 'SCSNMP2-CP-hatqz', 'SCSNMP2-CP-atqz',
            'SCSNMP2-CP-a5z', 'SCSNMP2-CP-aq5z'])
        self.plot_bars([None, None, None, None,
            'SCSMIMP2-CP-tz', 'SCSMIMP2-CP-matz', 'SCSMIMP2-CP-jatz', 'SCSMIMP2-CP-hatz', 'SCSMIMP2-CP-atz',
            'SCSMIMP2-CP-dtz', 'SCSMIMP2-CP-jadtz', 'SCSMIMP2-CP-hadtz', 'SCSMIMP2-CP-adtz',
            'SCSMIMP2-CP-qz', 'SCSMIMP2-CP-aaqz', 'SCSMIMP2-CP-maqz', 'SCSMIMP2-CP-jaqz', 'SCSMIMP2-CP-haqz', 'SCSMIMP2-CP-aqz',
            'SCSMIMP2-CP-tqz', 'SCSMIMP2-CP-matqz', 'SCSMIMP2-CP-jatqz', 'SCSMIMP2-CP-hatqz', 'SCSMIMP2-CP-atqz',
            None, None])
        self.plot_bars(['DWMP2-CP-dz', 'DWMP2-CP-jadz', 'DWMP2-CP-hadz', 'DWMP2-CP-adz',
            'DWMP2-CP-tz', 'DWMP2-CP-matz', 'DWMP2-CP-jatz', 'DWMP2-CP-hatz', 'DWMP2-CP-atz',
            'DWMP2-CP-dtz', 'DWMP2-CP-jadtz', 'DWMP2-CP-hadtz', 'DWMP2-CP-adtz',
            'DWMP2-CP-qz', 'DWMP2-CP-aaqz', 'DWMP2-CP-maqz', 'DWMP2-CP-jaqz', 'DWMP2-CP-haqz', 'DWMP2-CP-aqz',
            'DWMP2-CP-tqz', 'DWMP2-CP-matqz', 'DWMP2-CP-jatqz', 'DWMP2-CP-hatqz', 'DWMP2-CP-atqz',
            'DWMP2-CP-a5z', 'DWMP2-CP-aq5z'])
        self.plot_bars(['MP2C-CP-dz', 'MP2C-CP-jadz', 'MP2C-CP-hadz', 'MP2C-CP-adz',
            'MP2C-CP-tz', 'MP2C-CP-matz', 'MP2C-CP-jatz', 'MP2C-CP-hatz', 'MP2C-CP-atz',
            'MP2C-CP-dtz', 'MP2C-CP-jadtz', 'MP2C-CP-hadtz', 'MP2C-CP-adtz',
            None, None, None, None, None, 'MP2C-CP-aqz',
            None, None, None, None, 'MP2C-CP-atqz',
            None, None])
        self.plot_bars(['MP2C-CP-atqzdz', 'MP2C-CP-atqzjadz', 'MP2C-CP-atqzhadz', 'MP2C-CP-atqzadz',
            'MP2C-CP-atqztz', 'MP2C-CP-atqzmatz', 'MP2C-CP-atqzjatz', 'MP2C-CP-atqzhatz', 'MP2C-CP-atqzatz',
            'MP2C-CP-atqzdtz', 'MP2C-CP-atqzjadtz', 'MP2C-CP-atqzhadtz', 'MP2C-CP-atqzadtz'])

        # Fig. bars (c)
        self.plot_bars(['MP2F12-CP-dz', 'MP2F12-CP-jadz', 'MP2F12-CP-hadz', 'MP2F12-CP-adz',
            'MP2F12-CP-tz', 'MP2F12-CP-matz', 'MP2F12-CP-jatz', 'MP2F12-CP-hatz', 'MP2F12-CP-atz',
            'MP2F12-CP-dtz', 'MP2F12-CP-jadtz', 'MP2F12-CP-hadtz', 'MP2F12-CP-adtz',
            'MP2F12-CP-aqz', 'MP2F12-CP-atqz'])
        self.plot_bars(['SCSMP2F12-CP-dz', 'SCSMP2F12-CP-jadz', 'SCSMP2F12-CP-hadz', 'SCSMP2F12-CP-adz',
            'SCSMP2F12-CP-tz', 'SCSMP2F12-CP-matz', 'SCSMP2F12-CP-jatz', 'SCSMP2F12-CP-hatz', 'SCSMP2F12-CP-atz',
            'SCSMP2F12-CP-dtz', 'SCSMP2F12-CP-jadtz', 'SCSMP2F12-CP-hadtz', 'SCSMP2F12-CP-adtz',
            'SCSMP2F12-CP-aqz', 'SCSMP2F12-CP-atqz'])
        self.plot_bars(['SCSNMP2F12-CP-dz', 'SCSNMP2F12-CP-jadz', 'SCSNMP2F12-CP-hadz', 'SCSNMP2F12-CP-adz',
            'SCSNMP2F12-CP-tz', 'SCSNMP2F12-CP-matz', 'SCSNMP2F12-CP-jatz', 'SCSNMP2F12-CP-hatz', 'SCSNMP2F12-CP-atz',
            'SCSNMP2F12-CP-dtz', 'SCSNMP2F12-CP-jadtz', 'SCSNMP2F12-CP-adtz', 'SCSNMP2F12-CP-adtz',
            'SCSNMP2F12-CP-aqz', 'SCSNMP2F12-CP-atqz'])
        self.plot_bars([None, None, None, None,
            'SCSMIMP2F12-CP-tz', 'SCSMIMP2F12-CP-matz', 'SCSMIMP2F12-CP-jatz', 'SCSMIMP2F12-CP-hatz', 'SCSMIMP2F12-CP-atz',
            'SCSMIMP2F12-CP-dtz', 'SCSMIMP2F12-CP-jadtz', 'SCSMIMP2F12-CP-hadtz', 'SCSMIMP2F12-CP-adtz',
            'SCSMIMP2F12-CP-aqz', 'SCSMIMP2F12-CP-atqz'])
        self.plot_bars(['DWMP2F12-CP-dz', 'DWMP2F12-CP-jadz', 'DWMP2F12-CP-hadz', 'DWMP2F12-CP-adz',
            'DWMP2F12-CP-tz', 'DWMP2F12-CP-matz', 'DWMP2F12-CP-jatz', 'DWMP2F12-CP-hatz', 'DWMP2F12-CP-atz',
            'DWMP2F12-CP-dtz', 'DWMP2F12-CP-jadtz', 'DWMP2F12-CP-hadtz', 'DWMP2F12-CP-adtz',
            'DWMP2F12-CP-aqz', 'DWMP2F12-CP-atqz'])
        self.plot_bars(['MP2CF12-CP-dz', 'MP2CF12-CP-jadz', 'MP2CF12-CP-hadz', 'MP2CF12-CP-adz',
            'MP2CF12-CP-tz', 'MP2CF12-CP-matz', 'MP2CF12-CP-jatz', 'MP2CF12-CP-hatz', 'MP2CF12-CP-atz',
            'MP2CF12-CP-dtz', 'MP2CF12-CP-jadtz', 'MP2CF12-CP-hadtz', 'MP2CF12-CP-adtz',
            'MP2CF12-CP-aqz', 'MP2CF12-CP-atqz'])
        self.plot_bars(['MP2CF12-CP-atqzdz', 'MP2CF12-CP-atqzjadz', 'MP2CF12-CP-atqzhadz', 'MP2CF12-CP-atqzadz',
            'MP2CF12-CP-atqztz', 'MP2CF12-CP-atqzmatz', 'MP2CF12-CP-atqzjatz', 'MP2CF12-CP-atqzhatz', 'MP2CF12-CP-atqzatz',
            'MP2CF12-CP-atqzdtz', 'MP2CF12-CP-atqzjadtz', 'MP2CF12-CP-atqzhadtz', 'MP2CF12-CP-atqzadtz'])

    def make_pt2_Figure_2(self):
        """Plot all the graphics needed for the diffuse augmented grey
        bars plot in Fig. 2 of PT2.

        """
        # Fig. bars (a)
        self.plot_bars(['MP2-CP-adz', 'MP2-CP-atz', 'MP2-CP-adtz',
            'MP2-CP-aqz', 'MP2-CP-atqz', 'MP2-CP-a5z', 'MP2-CP-aq5z'])
        self.plot_bars(['SCSMP2-CP-adz', 'SCSMP2-CP-atz',
            'SCSMP2-CP-adtz', 'SCSMP2-CP-aqz', 'SCSMP2-CP-atqz',
            'SCSMP2-CP-a5z', 'SCSMP2-CP-aq5z'])
        self.plot_bars(['SCSNMP2-CP-adz', 'SCSNMP2-CP-atz',
            'SCSNMP2-CP-adtz', 'SCSNMP2-CP-aqz', 'SCSNMP2-CP-atqz',
            'SCSNMP2-CP-a5z', 'SCSNMP2-CP-aq5z'])
        self.plot_bars(['SCSMIMP2-CP-atz', 'SCSMIMP2-CP-atz',
            'SCSMIMP2-CP-adtz', 'SCSMIMP2-CP-aqz', 'SCSMIMP2-CP-atqz'])
        self.plot_bars(['SCSMIMP2-CP-tz', 'SCSMIMP2-CP-tz',
            'SCSMIMP2-CP-dtz', 'SCSMIMP2-CP-qz', 'SCSMIMP2-CP-tqz'])
        self.plot_bars(['DWMP2-CP-adz', 'DWMP2-CP-atz', 'DWMP2-CP-adtz',
            'DWMP2-CP-aqz', 'DWMP2-CP-atqz', 'DWMP2-CP-a5z', 'DWMP2-CP-aq5z'])
        self.plot_bars(['MP2C-CP-adz', 'MP2C-CP-adtzadz',
            'MP2C-CP-atqzadz', 'MP2C-CP-aq5zadz', 'MP2C-CP-atz',
            'MP2C-CP-atqzatz', 'MP2C-CP-aq5zatz', 'MP2C-CP-adtz',
            'MP2C-CP-atqzadtz', 'MP2C-CP-aqz', 'MP2C-CP-atqz'])

        # Fig. bars (b)
        self.plot_bars(['MP3-CP-adz', 'MP3-CP-adtzadz', 'MP3-CP-atqzadz',
            'MP3-CP-atz', 'MP3-CP-atqzatz', 'MP3-CP-adtz', 'MP3-CP-atqzadtz'])
        self.plot_bars(['MP25-CP-adz', 'MP25-CP-adtzadz', 'MP25-CP-atqzadz',
            'MP25-CP-atz', 'MP25-CP-atqzatz', 'MP25-CP-adtz', 'MP25-CP-atqzadtz'])
        self.plot_bars(['CCSD-CP-adz', 'CCSD-CP-adtzadz', 'CCSD-CP-atqzadz',
            'CCSD-CP-atz', 'CCSD-CP-atqzatz', 'CCSD-CP-adtz', 'CCSD-CP-atqzadtz'])
        self.plot_bars(['SCSCCSD-CP-adz', 'SCSCCSD-CP-adtzadz',
            'SCSCCSD-CP-atqzadz', 'SCSCCSD-CP-atz', 'SCSCCSD-CP-atqzatz',
            'SCSCCSD-CP-adtz', 'SCSCCSD-CP-atqzadtz'])
        self.plot_bars(['SCSMICCSD-CP-adz', 'SCSMICCSD-CP-adtzadz',
            'SCSMICCSD-CP-atqzadz', 'SCSMICCSD-CP-atz', 'SCSMICCSD-CP-atqzatz',
            'SCSMICCSD-CP-adtz', 'SCSMICCSD-CP-atqzadtz'])
        self.plot_bars(['CCSDT-CP-adz', 'CCSDT-CP-adtzadz',
            'CCSDT-CP-atqzadz', 'CCSDT-CP-atz', 'CCSDT-CP-atqzatz',
            'CCSDT-CP-adtz', 'CCSDT-CP-atqzadtz'])

        # Fig. bars (c)
        self.plot_bars(['MP2F12-CP-adz', 'MP2F12-CP-atz', 'MP2F12-CP-adtz',
            'MP2F12-CP-aqz', 'MP2F12-CP-atqz'])
        self.plot_bars(['SCSMP2F12-CP-adz', 'SCSMP2F12-CP-atz',
            'SCSMP2F12-CP-adtz', 'SCSMP2F12-CP-aqz', 'SCSMP2F12-CP-atqz'])
        self.plot_bars(['SCSNMP2F12-CP-adz', 'SCSNMP2F12-CP-atz',
            'SCSNMP2F12-CP-adtz', 'SCSNMP2F12-CP-aqz',
            'SCSNMP2F12-CP-atqz'])
        self.plot_bars(['SCSMIMP2F12-CP-atz', 'SCSMIMP2F12-CP-atz',
            'SCSMIMP2F12-CP-adtz', 'SCSMIMP2F12-CP-aqz',
            'SCSMIMP2F12-CP-atqz'])
        self.plot_bars(['SCSMIMP2F12-CP-tz', 'SCSMIMP2F12-CP-tz', 'SCSMIMP2F12-CP-dtz'])
        self.plot_bars(['DWMP2F12-CP-adz', 'DWMP2F12-CP-atz',
            'DWMP2F12-CP-adtz', 'DWMP2F12-CP-aqz', 'DWMP2F12-CP-atqz'])
        self.plot_bars(['MP2CF12-CP-adz', 'MP2CF12-CP-adtzadz',
            'MP2CF12-CP-atqzadz', 'MP2CF12-CP-atz', 'MP2CF12-CP-atqzatz',
            'MP2CF12-CP-adtz', 'MP2CF12-CP-atqzadtz', 'MP2CF12-CP-aqz',
            'MP2CF12-CP-atqz'])

        # Fig. bars (d)
        self.plot_bars(['CCSDAF12-CP-adz', 'CCSDAF12-CP-adtzadz', 'CCSDAF12-CP-atqzadz'])
        self.plot_bars(['CCSDBF12-CP-adz', 'CCSDBF12-CP-adtzadz', 'CCSDBF12-CP-atqzadz'])
        self.plot_bars(['SCSCCSDAF12-CP-adz', 'SCSCCSDAF12-CP-adtzadz', 'SCSCCSDAF12-CP-atqzadz'])
        self.plot_bars(['SCSCCSDBF12-CP-adz', 'SCSCCSDBF12-CP-adtzadz', 'SCSCCSDBF12-CP-atqzadz'])
        self.plot_bars(['SCMICCSDAF12-CP-adz', 'SCMICCSDAF12-CP-adtzadz', 'SCMICCSDAF12-CP-atqzadz'])
        self.plot_bars(['SCMICCSDBF12-CP-adz', 'SCMICCSDBF12-CP-adtzadz', 'SCMICCSDBF12-CP-atqzadz'])
        self.plot_bars(['CCSDTAF12-CP-adz', 'CCSDTAF12-CP-adtzadz', 'CCSDTAF12-CP-atqzadz'])
        self.plot_bars(['CCSDTBF12-CP-adz', 'CCSDTBF12-CP-adtzadz', 'CCSDTBF12-CP-atqzadz'])
        self.plot_bars(['DWCCSDTF12-CP-adz', 'DWCCSDTF12-CP-adtzadz', 'DWCCSDTF12-CP-atqzadz'])

    def table_generic(self, mtd, bas, columnplan, rowplan=['bas', 'mtd'],
        dbse=['DB4'], opt=['CP'], err=['mae'], sset=['tt'],
        benchmark='default', failoninc=True,
        landscape=False, standalone=True, subjoin=True,
        plotpath='', theme='', filename=None):
        """Prepares dictionary of errors for all combinations of *mtd*, *opt*,
        *bas* with respect to model chemistry *benchmark*, mindful of *failoninc*.
        Once error dictionary is ready, it and all other arguments are passed
        along to textables.table_generic.

        """
        # gather list of model chemistries for table
        mcs = ['-'.join(prod) for prod in itertools.product(mtd, opt, bas)]

        # compute errors
        serrors = {}
        for mc in mcs:
            serrors[mc] = {}
            for ss in self.sset.keys():
                errblock = self.compute_statistics(mc, benchmark=benchmark, sset=ss,
                    failoninc=failoninc, verbose=False, returnindiv=False)
                serrors[mc][ss] = {}
                serrors[mc][ss]['S22'] = None if errblock['S22'] is None else format_errors(errblock['S22'], mode=3)
                serrors[mc][ss]['NBC1'] = None if errblock['NBC1'] is None else format_errors(errblock['NBC1'], mode=3)
                serrors[mc][ss]['HBC1'] = None if errblock['HBC1'] is None else format_errors(errblock['HBC1'], mode=3)
                serrors[mc][ss]['HSG'] = None if errblock['HSG'] is None else format_errors(errblock['HSG'], mode=3)
                serrors[mc][ss]['DB4'] = format_errors(errblock['DB4'], mode=3)

        textables.table_generic(dbse=dbse, serrors=serrors,
            mtd=mtd, bas=bas, columnplan=columnplan, rowplan=rowplan,
            opt=opt, err=err, sset=sset,
            landscape=landscape, standalone=standalone, subjoin=subjoin,
            plotpath=plotpath, theme=theme, filename=filename)

    def table_merge_abbr(self, mtd, bas, opt=['CP'], err=['mae'], benchmark='default', failoninc=True, plotpath='analysis/flats/mplflat_', theme='smmerge'):
        """Specialization of table_generic into table with minimal statistics
        (three S22 and three overall) plus embedded slat diagram as suitable
        for main paper. A single table is formed in sections by *bas* with
        lines *mtd* within each section.

        """
        rowplan = ['bas', 'mtd']
        columnplan = [
            ['l', r"""Method \& Basis Set""", '', textables.label, {}],
            ['d', r'S22', 'HB', textables.val, {'ss': 'hb', 'db': 'S22'}],
            ['d', r'S22', 'MX/DD', textables.val, {'ss': 'mxdd', 'db': 'S22'}],
            ['d', r'S22', 'TT', textables.val, {'ss': 'tt', 'db': 'S22'}],
            ['d', r'Overall', 'HB', textables.val, {'ss': 'hb', 'db': 'DB4'}],
            ['d', r'Overall', 'MX/DD', textables.val, {'ss': 'mxdd', 'db': 'DB4'}],
            ['d', r'Overall', 'TT', textables.val, {'ss': 'tt', 'db': 'DB4'}],
            ['l', r"""Error Distribution\footnotemark[1]""", r"""\includegraphics[width=6.67cm,height=3.5mm]{%s%s.pdf}""" % (plotpath, 'blank'), textables.graphics, {}],
            ['d', r'Time', '', textables.val, {'ss': 'tt-5min', 'db': 'NBC1'}]]
            # TODO Time column not right at all

        self.table_generic(mtd=mtd, bas=bas, columnplan=columnplan, rowplan=rowplan,
            opt=opt, err=err,
            benchmark=benchmark, failoninc=failoninc,
            landscape=False, standalone=True, subjoin=True,
            plotpath=plotpath, theme=theme, filename=None)
        # TODO: not handled: filename, TODO switch standalone

    def table_merge_suppmat(self, mtd, bas, opt=['CP'], err=['mae'], benchmark='default', failoninc=True, plotpath='analysis/flats/mplflat_', theme='lgmerge'):
        """Specialization of table_generic into table with as many statistics
        as will fit (mostly fullcurve and a few 5min) plus embedded slat
        diagram as suitable for supplementary material. Multiple tables are
        formed, one for each in *bas* with lines *mtd* within each table.

        """
        rowplan = ['bas', 'mtd']
        columnplan = [
            ['l', r"""Method \& Basis Set""", '', textables.label, {}],
            ['d', 'S22', 'HB', textables.val, {'ss': 'hb', 'db': 'S22'}],
            ['d', 'S22', 'MX', textables.val, {'ss': 'mx', 'db': 'S22'}],
            ['d', 'S22', 'DD', textables.val, {'ss': 'dd', 'db': 'S22'}],
            ['d', 'S22', 'TT', textables.val, {'ss': 'tt', 'db': 'S22'}],
            ['d', 'NBC10', 'MX', textables.val, {'ss': 'mx', 'db': 'NBC1'}],
            ['d', 'NBC10', 'DD', textables.val, {'ss': 'dd', 'db': 'NBC1'}],
            ['d', 'NBC10', 'TT', textables.val, {'ss': 'tt', 'db': 'NBC1'}],
            ['d', 'HBC6', 'HB/TT', textables.val, {'ss': 'tt', 'db': 'HBC1'}],
            ['d', 'HSG', 'HB', textables.val, {'ss': 'hb', 'db': 'HSG'}],
            ['d', 'HSG', 'MX', textables.val, {'ss': 'mx', 'db': 'HSG'}],
            ['d', 'HSG', 'DD', textables.val, {'ss': 'dd', 'db': 'HSG'}],
            ['d', 'HSG', 'TT', textables.val, {'ss': 'tt', 'db': 'HSG'}],
            ['d', 'Avg', 'TT ', textables.val, {'ss': 'tt', 'db': 'DB4'}],
            ['l', r"""Error Distribution\footnotemark[1]""", r"""\includegraphics[width=6.67cm,height=3.5mm]{%s%s.pdf}""" % (plotpath, 'blank'), textables.graphics, {}],
            ['d', 'NBC10', r"""TT\footnotemark[2]""", textables.val, {'ss': 'tt-5min', 'db': 'NBC1'}],
            ['d', 'HBC6', r"""TT\footnotemark[2] """, textables.val, {'ss': 'tt-5min', 'db': 'HBC1'}],
            ['d', 'Avg', r"""TT\footnotemark[2]""", textables.val, {'ss': 'tt-5min', 'db': 'DB4'}]]

        self.table_generic(mtd=mtd, bas=bas, columnplan=columnplan, rowplan=rowplan,
            opt=opt, err=err,
            benchmark=benchmark, failoninc=failoninc,
            landscape=False, standalone=True, subjoin=False,
            plotpath=plotpath, theme=theme, filename=None)
        # TODO: not handled: filename, TODO switch standalone


class ThreeDatabases(object):
    """

    """
    def __init__(self):
        sys.path.append(os.path.dirname(__file__) + '/../databases')
        # S22 database
        self.s22 = Database('S22')
        # A24 database
        self.a24 = Database('A24')
        # HSG database
        self.hsg = Database('HSG')
        # subset assembly pattern
        self.sset = collections.OrderedDict()
        # assembly pattern for transspecies modelchems
        self.mc = {}

        # load up data and definitions
        self.load_pt2()
        self.load_dilabio()
        self.load_subsets()
        self.define_supersubsets()
        self.define_supermodelchems()

    def load_pt2(self):
        """Load qcdb.ReactionDatum results from standard location.

        """
        self.s22.load_pt2()
        self.hsg.load_pt2()

    def load_dilabio(self):
        """Load qcdb.ReactionDatum results from standard location.

        """
        self.a24.load_dilabio()
        self.a24.load_f12dilabio()

    def load_subsets(self):
        """Load subsets from standard generators.

        """
        self.s22.load_subsets()
        self.hsg.load_subsets()
        self.a24.load_subsets()

    def define_supersubsets(self):
        """

        """
        self.sset['tt'] = ['default', 'default', 'default']
        self.sset['hb'] = ['hb', 'hb', 'hb']
        self.sset['mx'] = ['mx', 'mx', 'mx']
        self.sset['dd'] = ['dd', 'dd', 'dd']
        self.sset['mxdd'] = ['mxdd', 'mxdd', 'mxdd']
        self.sset['pp'] = ['mxddpp', 'mxddpp', 'mxddpp']
        self.sset['np'] = ['mxddnp', 'mxddnp', 'mxddnp']
        self.sset['tt-5min'] = ['default', 'default', 'default']
        self.sset['hb-5min'] = ['hb', 'hb', 'hb']
        self.sset['mx-5min'] = ['mx', 'mx', 'mx']
        self.sset['dd-5min'] = ['dd', 'dd', 'dd']
        self.sset['mxdd-5min'] = ['mxdd', 'mxdd', 'mxdd']
        self.sset['pp-5min'] = ['mxddpp', 'mxddpp', 'mxddpp']
        self.sset['np-5min'] = ['mxddnp', 'mxddnp', 'mxddnp']
        self.sset['weak'] = ['weak', 'weak', 'weak']
        self.sset['weak_hb'] = ['weak_hb', None, 'weak_hb']
        self.sset['weak_mx'] = ['weak_mx', 'weak_mx', 'weak_mx']
        self.sset['weak_dd'] = ['weak_dd', 'weak_dd', 'weak_dd']

    def define_supermodelchems(self):
        """

        """
        self.mc['CCSD-CP-adz'] = ['CCSD-CP-adz', 'CCSD-CP-hadz', 'CCSD-CP-adz']
        self.mc['CCSD-CP-atz'] = ['CCSD-CP-atz', 'CCSD-CP-hatz', 'CCSD-CP-atz']
        self.mc['CCSD-CP-adtz'] = ['CCSD-CP-adtz', 'CCSD-CP-hadtz', 'CCSD-CP-adtz']
        self.mc['CCSD-CP-adtzadz'] = ['CCSD-CP-adtzadz', 'CCSD-CP-adtzhadz', 'CCSD-CP-adtzadz']
        self.mc['CCSD-CP-atzadz'] = ['CCSD-CP-atzadz', 'CCSD-CP-atzhadz', 'CCSD-CP-atzadz']
        self.mc['CCSD-CP-atqzadz'] = ['CCSD-CP-atqzadz', 'CCSD-CP-atqzhadz', 'CCSD-CP-atqzadz']
        self.mc['CCSD-CP-atzadtz'] = ['CCSD-CP-atzadtz', 'CCSD-CP-atzhadtz', 'CCSD-CP-atzadtz']
        self.mc['CCSD-CP-atqzadtz'] = ['CCSD-CP-atqzadtz', 'CCSD-CP-atqzhadtz', 'CCSD-CP-atqzadtz']
        self.mc['CCSD-CP-atqzatz'] = ['CCSD-CP-atqzatz', 'CCSD-CP-atqzhatz', 'CCSD-CP-atqzatz']

        #self.mc['SCSCCSD-CP-adz'] = ['SCSCCSD-CP-adz', 'SCSCCSD-CP-hadz', 'SCSCCSD-CP-adz', 'SCSCCSD-CP-hadz']
        #self.mc['SCSCCSD-CP-atz'] = ['SCSCCSD-CP-atz', 'SCSCCSD-CP-hatz', 'SCSCCSD-CP-atz', 'SCSCCSD-CP-hatz']
        #self.mc['SCSCCSD-CP-adtz'] = ['SCSCCSD-CP-adtz', 'SCSCCSD-CP-hadtz', 'SCSCCSD-CP-adtz', 'SCSCCSD-CP-hadtz']
        #self.mc['SCSCCSD-CP-adtzadz'] = ['SCSCCSD-CP-adtzadz', 'SCSCCSD-CP-adtzhadz', 'SCSCCSD-CP-adtzadz', 'SCSCCSD-CP-adtzhadz']
        #self.mc['SCSCCSD-CP-atzadz'] = ['SCSCCSD-CP-atzadz', 'SCSCCSD-CP-atzhadz', 'SCSCCSD-CP-atzadz', 'SCSCCSD-CP-atzhadz']
        #self.mc['SCSCCSD-CP-atqzadz'] = ['SCSCCSD-CP-atqzadz', 'SCSCCSD-CP-atqzhadz', 'SCSCCSD-CP-atqzadz', 'SCSCCSD-CP-atqzhadz']
        #self.mc['SCSCCSD-CP-atzadtz'] = ['SCSCCSD-CP-atzadtz', 'SCSCCSD-CP-atzhadtz', 'SCSCCSD-CP-atzadtz', 'SCSCCSD-CP-atzhadtz']
        #self.mc['SCSCCSD-CP-atqzadtz'] = ['SCSCCSD-CP-atqzadtz', 'SCSCCSD-CP-atqzhadtz', 'SCSCCSD-CP-atqzadtz', 'SCSCCSD-CP-atqzhadtz']
        #self.mc['SCSCCSD-CP-atqzatz'] = ['SCSCCSD-CP-atqzatz', 'SCSCCSD-CP-atqzhatz', 'SCSCCSD-CP-atqzatz', 'SCSCCSD-CP-atqzhatz']

        #self.mc['SCSMICCSD-CP-adz'] = ['SCSMICCSD-CP-adz', 'SCSMICCSD-CP-hadz', 'SCSMICCSD-CP-adz', 'SCSMICCSD-CP-hadz']
        #self.mc['SCSMICCSD-CP-atz'] = ['SCSMICCSD-CP-atz', 'SCSMICCSD-CP-hatz', 'SCSMICCSD-CP-atz', 'SCSMICCSD-CP-hatz']
        #self.mc['SCSMICCSD-CP-adtz'] = ['SCSMICCSD-CP-adtz', 'SCSMICCSD-CP-hadtz', 'SCSMICCSD-CP-adtz', 'SCSMICCSD-CP-hadtz']
        #self.mc['SCSMICCSD-CP-adtzadz'] = ['SCSMICCSD-CP-adtzadz', 'SCSMICCSD-CP-adtzhadz', 'SCSMICCSD-CP-adtzadz', 'SCSMICCSD-CP-adtzhadz']
        #self.mc['SCSMICCSD-CP-atzadz'] = ['SCSMICCSD-CP-atzadz', 'SCSMICCSD-CP-atzhadz', 'SCSMICCSD-CP-atzadz', 'SCSMICCSD-CP-atzhadz']
        #self.mc['SCSMICCSD-CP-atqzadz'] = ['SCSMICCSD-CP-atqzadz', 'SCSMICCSD-CP-atqzhadz', 'SCSMICCSD-CP-atqzadz', 'SCSMICCSD-CP-atqzhadz']
        #self.mc['SCSMICCSD-CP-atzadtz'] = ['SCSMICCSD-CP-atzadtz', 'SCSMICCSD-CP-atzhadtz', 'SCSMICCSD-CP-atzadtz', 'SCSMICCSD-CP-atzhadtz']
        #self.mc['SCSMICCSD-CP-atqzadtz'] = ['SCSMICCSD-CP-atqzadtz', 'SCSMICCSD-CP-atqzhadtz', 'SCSMICCSD-CP-atqzadtz', 'SCSMICCSD-CP-atqzhadtz']
        #self.mc['SCSMICCSD-CP-atqzatz'] = ['SCSMICCSD-CP-atqzatz', 'SCSMICCSD-CP-atqzhatz', 'SCSMICCSD-CP-atqzatz', 'SCSMICCSD-CP-atqzhatz']

        self.mc['CCSDT-CP-adz'] = ['CCSDT-CP-adz', 'CCSDT-CP-hadz', 'CCSDT-CP-adz']
        self.mc['CCSDT-CP-atz'] = ['CCSDT-CP-atz', 'CCSDT-CP-hatz', 'CCSDT-CP-atz']
        self.mc['CCSDT-CP-adtz'] = ['CCSDT-CP-adtz', 'CCSDT-CP-hadtz', 'CCSDT-CP-adtz']
        self.mc['CCSDT-CP-adtzadz'] = ['CCSDT-CP-adtzadz', 'CCSDT-CP-adtzhadz', 'CCSDT-CP-adtzadz']
        self.mc['CCSDT-CP-atzadz'] = ['CCSDT-CP-atzadz', 'CCSDT-CP-atzhadz', 'CCSDT-CP-atzadz']
        self.mc['CCSDT-CP-atqzadz'] = ['CCSDT-CP-atqzadz', 'CCSDT-CP-atqzhadz', 'CCSDT-CP-atqzadz']
        self.mc['CCSDT-CP-atzadtz'] = ['CCSDT-CP-atzadtz', 'CCSDT-CP-atzhadtz', 'CCSDT-CP-atzadtz']
        self.mc['CCSDT-CP-atqzadtz'] = ['CCSDT-CP-atqzadtz', 'CCSDT-CP-atqzhadtz', 'CCSDT-CP-atqzadtz']
        self.mc['CCSDT-CP-atqzatz'] = ['CCSDT-CP-atqzatz', 'CCSDT-CP-atqzhatz', 'CCSDT-CP-atqzatz']

    def analyze_modelchems(self, modelchem, benchmark='default', failoninc=True, verbose=False):
        """Print out nicely formatted summary statistics for every model
        chemistry in array *modelchem* versus *benchmark* for every
        registered subset.

        """
        # compute errors
        errors = {}
        for mc in modelchem:
            errors[mc] = {}
            for ss in self.sset.keys():
                errors[mc][ss] = self.compute_statistics(mc, benchmark=benchmark, sset=ss,
                    failoninc=failoninc, verbose=verbose, returnindiv=False)
        # present errors
        pre, suf, mid = string_contrast(modelchem)
        print """\n  ==> %s %s[]%s Errors <==""" % ('DB3', pre, suf)
        print """%20s    %42s%42s%42s%42s""" % \
            ('', '=> DB3 <=', '=> S22 <=', '=> HSG <=', '=> A24 <=')
        print """%20s        %5s  %4s   %6s %6s    %6s""" % \
            ('', 'ME', 'STDE', 'MAE', 'MA%E', 'MA%BE')
        for ss in self.sset.keys():
            print """   => %s <= """ % (ss)
            for mc in modelchem:
                tmpdb = errors[mc][ss]
                print """%20s    %42s%42s%42s%42s""" % (mid[modelchem.index(mc)],
                    format_errors(errors[mc][ss]['DB3']),
                    '' if tmpdb['S22'] is None else format_errors(tmpdb['S22']),
                    '' if tmpdb['HSG'] is None else format_errors(tmpdb['HSG']),
                    '' if tmpdb['A24'] is None else format_errors(tmpdb['A24']))

    def plot_bars(self, modelchem, benchmark='default', sset=['tt-5min', 'hb-5min', 'mx-5min', 'dd-5min'], failoninc=True, verbose=False):
        """Prepares 'grey bars' diagram for each model chemistry in array
        *modelchem* versus *benchmark* over all four databases. A wide bar
        is plotted with three smaller bars, corresponding to the 'mae'
        summary statistic of the four subsets in *sset*. Prepares bars
        diagram instructions and either executes them if matplotlib
        available (Canopy) or prints them.

        """
        # compute errors
        errors = {}
        for mc in modelchem:
            if mc is not None:
                errors[mc] = {}
                for ss in sset:
                    errors[mc][ss] = self.compute_statistics(mc, benchmark=benchmark, sset=ss,
                        failoninc=failoninc, verbose=verbose, returnindiv=False)
        # repackage
        pre, suf, mid = string_contrast(modelchem)
        #dbdat = [{'mc': mid[modelchem.index(mc)], 'data': [errors[mc][ss]['DB4']['mae'] for ss in sset]} for mc in modelchem]
        dbdat = []
        for mc in modelchem:
            if mc is None:
                dbdat.append(None)
            else:
                dbdat.append({'mc': mid[modelchem.index(mc)], 'data': [errors[mc][ss]['DB3']['mae'] for ss in sset]})
        title = '3DB ' + pre + '[]' + suf
        # generate matplotlib instructions and call or print
        try:
            import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print """mpl.bar(%s,\n    title='%s')\n\n""" % (dbdat, title)
        else:
            # if running from Canopy, call mpl directly
            mpl.bar(dbdat, title=title)

    def compute_statistics(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False, returnindiv=False):
        """Computes summary statistics and, if *returnindiv* True,
        individual errors for single model chemistry *modelchem* versus
        *benchmark* over subset *sset* over all four databases.
        Particularly, imposes cross-database definitions for sset and
        modelchem.

        """
        pre, suf, mid = string_contrast(modelchem)
        if modelchem in self.mc.keys():
            lmc = self.mc[modelchem]
        else:
            lmc = [modelchem, modelchem, modelchem, modelchem]
        if sset in self.sset.keys():
            lss = self.sset[sset]
        else:
            lss = [sset, sset, sset, sset]

        errors = collections.OrderedDict()
        indiv = collections.OrderedDict()
        errors['S22'], indiv['S22'] = (None, None) if lss[0] is None else self.s22.compute_statistics(lmc[0], sset=lss[0],
            benchmark=benchmark, failoninc=failoninc, verbose=verbose, returnindiv=True)
        errors['HSG'], indiv['HSG'] = (None, None) if lss[1] is None else self.hsg.compute_statistics(lmc[1], sset=lss[1],
            benchmark=benchmark, failoninc=failoninc, verbose=verbose, returnindiv=True)
        errors['A24'], indiv['A24'] = (None, None) if lss[2] is None else self.a24.compute_statistics(lmc[2], sset=lss[2],
            benchmark=benchmark, failoninc=failoninc, verbose=verbose, returnindiv=True)

        args = []
        if lss[0] is not None:
            args.append(errors['S22'])
        if lss[1] is not None:
            args.append(errors['HSG'])
        if lss[2] is not None:
            args.append(errors['A24'])
        errors['DB3'] = average_errors(*args)

        if returnindiv:
            return errors, indiv
        else:
            return errors

    def plot_modelchems(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False, color='sapt', xlimit=4.0):
        """Computes individual errors and summary statistics for each
        model chemistry in array *modelchem* versus *benchmark* over
        subset *sset* over all four databases. Thread *color* can be 'rgb'
        for old coloring, a color name or 'sapt' for spectrum coloring.
        Prepares thread diagram instructions and either executes them if
        matplotlib available (Canopy) or prints them.

        """
        # compute errors
        errors = {}
        indiv = {}
        for mc in modelchem:
            errors[mc], indiv[mc] = self.compute_statistics(mc, benchmark=benchmark, sset=sset,
                failoninc=failoninc, verbose=verbose, returnindiv=True)
        # repackage
        dbdat = []
        for db in indiv[modelchem[0]].keys():
            if indiv[modelchem[0]][db] is not None:
                for rxn, orxn in indiv[modelchem[0]][db].items():
                    if db == 'S22':
                        lcolor = self.s22.hrxn[rxn].color
                    elif db == 'HSG':
                        lcolor = self.hsg.hrxn[rxn].color
                    elif db == 'A24':
                        lcolor = self.a24.hrxn[rxn].color
                    dbdat.append({'sys': str(rxn), 'color': lcolor,
                        'data': [indiv[mc][db][rxn][0] for mc in modelchem]})
        pre, suf, mid = string_contrast(modelchem)
        title = '3DB-' + sset + ' ' + pre + '[]' + suf
        mae = [errors[mc]['DB3']['mae'] for mc in modelchem]
        mapbe = [100 * errors[mc]['DB3']['mapbe'] for mc in modelchem]
        # generate matplotlib instructions and call or print
        try:
            import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print """mpl.thread(%s,\n    color='%s',\n    title='%s',\n    labels=%s,\n    mae=%s,\n    mape=%s\n    xlimit=%s)\n\n""" % \
                (dbdat, color, title, mid, mae, mapbe, str(xlimit))
        else:
            # if running from Canopy, call mpl directly
            mpl.thread(dbdat, color=color, title=title, labels=mid, mae=mae, mape=mapbe, xlimit=xlimit)

    def plot_flat(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False, color='sapt', xlimit=4.0, view=True):
        """Computes individual errors and summary statistics for single
        model chemistry *modelchem* versus *benchmark* over
        subset *sset* over all four databases. Thread *color* can be 'rgb'
        for old coloring, a color name or 'sapt' for spectrum coloring.
        Prepares flat diagram instructions and either executes them if
        matplotlib available (Canopy) or prints them.

        """
        # compute errors
        errors = {}
        indiv = {}
        mc = modelchem
        errors[mc], indiv[mc] = self.compute_statistics(mc, benchmark=benchmark, sset=sset,
            failoninc=failoninc, verbose=verbose, returnindiv=True)
        # repackage
        dbdat = []
        for db in indiv[mc].keys():
            if indiv[mc][db] is not None:
                for rxn, orxn in indiv[mc][db].items():
                    if db == 'S22':
                        lcolor = self.s22.hrxn[rxn].color
                    elif db == 'HSG':
                        lcolor = self.hsg.hrxn[rxn].color
                    elif db == 'A24':
                        lcolor = self.a24.hrxn[rxn].color
                    dbdat.append({'sys': str(rxn), 'color': lcolor,
                        'data': [indiv[mc][db][rxn][0]]})
        pre, suf, mid = string_contrast(mc)
        title = '3DB-' + sset + ' ' + pre + '[]' + suf
        mae = errors[mc]['DB3']['mae']
        mapbe = 100 * errors[mc]['DB3']['mapbe']
        mapbe = None
        # generate matplotlib instructions and call or print
        try:
            import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print """mpl.flat(%s,\n    color='%s',\n    title='%s',\n    mae=%s,\n    mape=%s,\n    xlimit=%s,\n    view=%s)\n\n""" % \
                (dbdat, color, mc, mae, mapbe, xlimit, view)
        else:
            # if running from Canopy, call mpl directly
            mpl.flat(dbdat, color=color, title=mc, mae=mae, mape=mapbe, xlimit=xlimit, view=view)

# print certain statistic for all 4 db and summary and indev sys if min or max
#asdf.analyze_modelchems(['CCSD-CP-adz', 'CCSD-CP-adtzadz', 'CCSD-CP-atqzadz', 'CCSD-CP-atz', 'CCSD-CP-atqzatz', 'CCSD-CP-adtz', 'CCSD-CP-atqzadtz'])
