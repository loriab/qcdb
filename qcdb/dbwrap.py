import os
import sys
import math
import collections
from exceptions import *
from modelchems import Method, BasisSet, methods, bases
import psiutil


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
    elif mode == 2:
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


def string_contrast(s):
    """From an array of strings, *s*, returns maximum common prefix
    string, maximum common suffix string, and array of middles.

    """
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

    middle = [mc[preidx:sufidx] for mc in s]
    prefix = short[:preidx]
    suffix = short[sufidx:]

    return prefix, suffix, middle


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

#        self.sset = {}
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
        statistics between *modelchem* and *benchmark* model
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

    def load_pt2(self, modname=None, funcname='load_pt2', pythonpath=None):
        """Loads qcdb.ReactionDatums from module *modname* function
        *funcname* (which default to self.dbse + '_pt2'
        and 'load_pt2').

        """
        if pythonpath is not None:
            sys.path.insert(0, pythonpath)
        else:
            sys.path.append(os.path.dirname(__file__) + '/../data')
        try:
            pt2mod = __import__(self.dbse + '_pt2') if modname is None else __import__(modname)
        except ImportError:
            if modname is None:
                print """Wavefunction/PT2 data unavailable for database %s.\n""" % (self.dbse)
                return
            else:
                print '\nPython module for database data %s failed to load\n\n' % (modname)
                print '\nSearch path that was tried:\n'
                print ", ".join(map(str, sys.path))
                raise ValidationError("Python module loading problem for database data " + str(modname))
        try:
            getattr(pt2mod, funcname)(self)
        except AttributeError:
            raise ValidationError("Python module missing function %s for loading data " % (str(funcname)))

        print """Database %s: Wavefunction results loaded""" % (self.dbse)

    def load_dhdft(self, modname=None, funcname='load_dhdft', pythonpath=None):
        """Loads qcdb.ReactionDatums from module *modname* function
        *funcname* (which default to self.dbse + '_dhdft'
        and 'load_dhdft').
 
        """
        if pythonpath is not None:
            sys.path.insert(0, pythonpath)
        else:
            sys.path.append(os.path.dirname(__file__) + '/../data')
        try:
            pt2mod = __import__(self.dbse + '_dhdft') if modname is None else __import__(modname)
        except ImportError:
            if modname is None:
                print """DH-DFT data unavailable for database %s.\n""" % (self.dbse)
                return
            else:
                print '\nPython module for database data %s failed to load\n\n' % (modname)
                print '\nSearch path that was tried:\n'
                print ", ".join(map(str, sys.path))
                raise ValidationError("Python module loading problem for database data " + str(modname))
        try:
            getattr(pt2mod, funcname)(self)
        except AttributeError:
            raise ValidationError("Python module missing function %s for loading data " % (str(funcname)))
 
        print """Database %s: DH-DFT results loaded""" % (self.dbse)

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

    def analyze_modelchem(self, modelchem, benchmark='default', failoninc=True, verbose=False):
        """Compute and print error statistics for *modelchem* versus
        *benchmark* for all available subsets and return dictionary of same.

        """
        # compute errors
        errors = collections.OrderedDict()
        for ss in self.sset.keys():
            errors[ss] = self.compute_statistics(modelchem, benchmark=benchmark, sset=ss, failoninc=failoninc, verbose=verbose)
        # present errors
        print """\n  ==> %s %s Errors <==""" % (self.dbse, modelchem)
        print """%20s        %5s  %4s   %6s %6s    %6s""" % \
            ('', 'ME', 'STDE', 'MAE', 'MA%E', 'MA%BE')
        for ss in errors.keys():
            if any(errors[ss].values()):
                print """%20s    %42s""" % (ss, format_errors(errors[ss]))
        # return errors
        return errors

    def analyze_modelchems(self, modelchem, benchmark='default', failoninc=True, verbose=False):
        """Compute and print error statistics for each model chemistry in
        array *modelchem* versus *benchmark* for all available subsets and
        return dictionary of same.

        """
        pre, suf, mid = string_contrast(modelchem)
        # compute errors
        errors = collections.OrderedDict()
        for ss in self.sset.keys():
            errors[ss] = collections.OrderedDict()
            for mc in modelchem:
                errors[ss][mc] = self.compute_statistics(mc, benchmark=benchmark, sset=ss, failoninc=failoninc, verbose=verbose)
        # present errors
        print """\n  ==> %s %s[]%s Errors <==""" % (self.dbse, pre, suf)
        print """%20s        %5s  %4s   %6s %6s    %6s""" % \
            ('', 'ME', 'STDE', 'MAE', 'MA%E', 'MA%BE')
        for ss in self.sset.keys():
            if any([any(errors[ss][mc].values()) for mc in modelchem]):
                print """   => %s <= """ % (ss)
                for mc in modelchem:
                    print """%20s    %42s""" % (mid[modelchem.index(mc)], format_errors(errors[ss][mc]))
        # return errors
        return errors

    def plot_modelchems(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False, color='sapt'):
        """Computes individual errors and summary statistics for *modelchem*
        versus *benchmark* over subset *sset*.  Thread *color* can be 'rgb' for
        old coloring, a color name or 'sapt' for spectrum coloring.  Prepares
        thread diagram instructions and either executes them if matplotlib
        available (Canopy) or prints them.

        """
        pre, suf, mid = string_contrast(modelchem)
        # compute errors
        errors = collections.OrderedDict()
        indiv = collections.OrderedDict()
        for mc in modelchem:
            errors[mc], indiv[mc] = self.compute_statistics(mc, benchmark=benchmark,
                sset=sset, failoninc=failoninc, verbose=verbose, returnindiv=True)
        # repackage
        dbdat = [{'dbse': self.dbse, 'sys': str(rxn), 'color': self.hrxn[rxn].color, 'data': [indiv[mc][rxn][0] for mc in modelchem]} for rxn in self.sset[sset].keys()]
        title = self.dbse + ' ' + pre + '[]' + suf
        mae = [errors[mc]['mae'] for mc in modelchem]
        mapbe = [100 * errors[mc]['mapbe'] for mc in modelchem]
        # generate matplotlib instructions and call or print
        try:
            import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print """mpl.thread(%s,\n    color='%s',\n    title='%s',\n    labels=%s,\n    mae=%s,\n    mape=%s)\n\n""" % \
                (dbdat, color, title, mid, mae, mapbe)
        else:
            # if running from Canopy, call mpl directly
            mpl.thread(dbdat, color=color, title=title, labels=mid, mae=mae, mape=mapbe)


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

        # load up data and definitions
        self.load_pt2()
        self.load_dhdft()
        self.load_subsets()
        self.define_supersubsets()

    def load_pt2(self):
        """Load qcdb.ReactionDatum results from standard location.

        """
        self.s22.load_pt2()
        self.nbc1.load_pt2()
        self.hbc1.load_pt2()
        self.hsg.load_pt2()

    def load_dhdft(self):
        """Load qcdb.ReactionDatum results from standard location.

        """
        self.s22.load_dhdft()
        self.nbc1.load_dhdft()
        self.hbc1.load_dhdft()
        self.hsg.load_dhdft()

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

    def analyze_modelchem(self, modelchem, benchmark='default', failoninc=True, verbose=False):
        """

        """
        errors = collections.OrderedDict()
        errors['S22'] = self.s22.analyze_modelchem(modelchem, benchmark=benchmark, failoninc=failoninc, verbose=verbose)
        errors['NBC1'] = self.nbc1.analyze_modelchem(modelchem, benchmark=benchmark, failoninc=failoninc, verbose=verbose)
        errors['HBC1'] = self.hbc1.analyze_modelchem(modelchem, benchmark=benchmark, failoninc=failoninc, verbose=verbose)
        errors['HSG'] = self.hsg.analyze_modelchem(modelchem, benchmark=benchmark, failoninc=failoninc, verbose=verbose)
        errors['DB4'] = collections.OrderedDict()

        for ss, arr in self.sset.items():
            args = []
            if arr[0] is not None:
                args.append(errors['S22'][arr[0]])
            if arr[1] is not None:
                args.append(errors['NBC1'][arr[1]])
            if arr[2] is not None:
                args.append(errors['HBC1'][arr[2]])
            if arr[3] is not None:
                args.append(errors['HSG'][arr[3]])
            errors['DB4'][ss] = average_errors(*args)

        print """\n  ==> %s %s Errors <==""" % ('DB4', modelchem)
        print """%20s    %42s%42s%42s%42s%42s""" % \
            ('', '=> DB4 <=', '=> S22 <=', '=> NBC1 <=', '=> HBC1 <=', '=> HSG <=')
        print """%20s        %5s  %4s   %6s %6s    %6s""" % \
            ('', 'ME', 'STDE', 'MAE', 'MA%E', 'MA%BE')
        for ss, arr in self.sset.items():
            print """%20s    %42s%42s%42s%42s%42s""" % (ss, format_errors(errors['DB4'][ss]),
                '' if arr[0] is None else format_errors(errors['S22'][arr[0]]),
                '' if arr[1] is None else format_errors(errors['NBC1'][arr[1]]),
                '' if arr[2] is None else format_errors(errors['HBC1'][arr[2]]),
                '' if arr[3] is None else format_errors(errors['HSG'][arr[3]]))

    def analyze_modelchems(self, modelchem, benchmark='default', failoninc=True, verbose=False):
        """

        """
        pre, suf, mid = string_contrast(modelchem)
        # compute errors
        errors = collections.OrderedDict()
        errors['S22'] = self.s22.analyze_modelchems(modelchem, benchmark=benchmark, failoninc=failoninc, verbose=verbose)
        errors['NBC1'] = self.nbc1.analyze_modelchems(modelchem, benchmark=benchmark, failoninc=failoninc, verbose=verbose)
        errors['HBC1'] = self.hbc1.analyze_modelchems(modelchem, benchmark=benchmark, failoninc=failoninc, verbose=verbose)
        errors['HSG'] = self.hsg.analyze_modelchems(modelchem, benchmark=benchmark, failoninc=failoninc, verbose=verbose)
        errors['DB4'] = collections.OrderedDict()
        for ss, arr in self.sset.items():
            errors['DB4'][ss] = collections.OrderedDict()
            for mc in modelchem:
                args = []
                if arr[0] is not None:
                    args.append(errors['S22'][arr[0]][mc])
                if arr[1] is not None:
                    args.append(errors['NBC1'][arr[1]][mc])
                if arr[2] is not None:
                    args.append(errors['HBC1'][arr[2]][mc])
                if arr[3] is not None:
                    args.append(errors['HSG'][arr[3]][mc])
                errors['DB4'][ss][mc] = average_errors(*args)
        # present errors
        print """\n  ==> %s %s[]%s Errors <==""" % ('DB4', pre, suf)
        print """%20s    %42s%42s%42s%42s%42s""" % \
            ('', '=> DB4 <=', '=> S22 <=', '=> NBC1 <=', '=> HBC1 <=', '=> HSG <=')
        print """%20s        %5s  %4s   %6s %6s    %6s""" % \
            ('', 'ME', 'STDE', 'MAE', 'MA%E', 'MA%BE')
        for ss, arr in self.sset.items():
            print """   => %s <= """ % (ss)
            for mc in modelchem:
                print """%20s    %42s%42s%42s%42s%42s""" % (mid[modelchem.index(mc)],
                    format_errors(errors['DB4'][ss][mc]),
                    '' if arr[0] is None else format_errors(errors['S22'][arr[0]][mc]),
                    '' if arr[1] is None else format_errors(errors['NBC1'][arr[1]][mc]),
                    '' if arr[2] is None else format_errors(errors['HBC1'][arr[2]][mc]),
                    '' if arr[3] is None else format_errors(errors['HSG'][arr[3]][mc]))

    def plot_modelchems(self, modelchem, benchmark='default', sset='default', failoninc=True, verbose=False, color='sapt'):
        """Computes individual errors and summary statistics for *modelchem*
        versus *benchmark* over subset *sset* over all four databases. Thread
        *color* can be 'rgb' for old coloring, a color name or 'sapt' for spectrum
        coloring.  Prepares thread diagram instructions and either executes them
        if matplotlib available (Canopy) or prints them.

        """
        pre, suf, mid = string_contrast(modelchem)
        # compute errors
        errors = collections.OrderedDict()
        indiv = collections.OrderedDict()
        errors = {}
        indiv = {}
        sset = sset.lower()
        indiv['S22'] = {} #collections.OrderedDict()
        indiv['NBC1'] = {} #collections.OrderedDict()
        indiv['HBC1'] = {} #collections.OrderedDict()
        indiv['HSG'] = {} #collections.OrderedDict()
        for mc in modelchem:
#            indiv[mc] = {}
            args = []
            if self.sset[sset][0] is not None:
                tmpe, indiv['S22'][mc] = self.s22.compute_statistics(mc, benchmark=benchmark,
                    sset=self.sset[sset][0], failoninc=failoninc, verbose=verbose, returnindiv=True)
                args.append(tmpe)
            if self.sset[sset][1] is not None:
                tmpe, indiv['NBC1'][mc] = self.nbc1.compute_statistics(mc, benchmark=benchmark,
                    sset=self.sset[sset][1], failoninc=failoninc, verbose=verbose, returnindiv=True)
                args.append(tmpe)
            if self.sset[sset][2] is not None:
                tmpe, indiv['HBC1'][mc] = self.hbc1.compute_statistics(mc, benchmark=benchmark,
                    sset=self.sset[sset][2], failoninc=failoninc, verbose=verbose, returnindiv=True)
                args.append(tmpe)
            if self.sset[sset][3] is not None:
                tmpe, indiv['HSG'][mc] = self.hsg.compute_statistics(mc, benchmark=benchmark,
                    sset=self.sset[sset][3], failoninc=failoninc, verbose=verbose, returnindiv=True)
                args.append(tmpe)
            errors[mc] = average_errors(*args)

        # repackage
        dbdat = []
        if self.sset[sset][0] is not None:
            for rxn in self.s22.sset[self.sset[sset][0]].keys():
                dbdat.append({'dbse': 'S22',
                              'sys': str(rxn),
                              'color': self.s22.hrxn[rxn].color,
                              'data': [indiv['S22'][mc][rxn][0] for mc in modelchem]}) 
        if self.sset[sset][1] is not None:
            for rxn in self.nbc1.sset[self.sset[sset][1]].keys():
                dbdat.append({'dbse': 'NBC1',
                              'sys': str(rxn),
                              'color': self.nbc1.hrxn[rxn].color,
                              'data': [indiv['NBC1'][mc][rxn][0] for mc in modelchem]}) 
        if self.sset[sset][2] is not None:
            for rxn in self.hbc1.sset[self.sset[sset][2]].keys():
                dbdat.append({'dbse': 'HBC1',
                              'sys': str(rxn),
                              'color': self.hbc1.hrxn[rxn].color,
                              'data': [indiv['HBC1'][mc][rxn][0] for mc in modelchem]}) 
        if self.sset[sset][3] is not None:
            for rxn in self.hsg.sset[self.sset[sset][3]].keys():
                dbdat.append({'dbse': 'HSG',
                              'sys': str(rxn),
                              'color': self.hsg.hrxn[rxn].color,
                              'data': [indiv['HSG'][mc][rxn][0] for mc in modelchem]}) 
#        dbdat = [{'dbse': self.dbse, 'sys': str(rxn), 'data': [indiv[mc][rxn][0] for mc in modelchem]} for rxn in self.sset[sset].keys()]
        title = '4DB-' + sset + ' ' + pre + '[]' + suf
        mae = [errors[mc]['mae'] for mc in modelchem]
        mapbe = [100 * errors[mc]['mapbe'] for mc in modelchem]
        # generate matplotlib instructions and call or print
        try:
            import mpl
            import matplotlib.pyplot as plt
        except ImportError:
            # if not running from Canopy, print line to execute from Canopy
            print """mpl.thread(%s,\n    color='%s',\n    title='%s',\n    labels=%s,\n    mae=%s,\n    mape=%s)\n\n""" % \
                (dbdat, color, title, mid, mae, mapbe)
        else:
            # if running from Canopy, call mpl directly
            mpl.thread(dbdat, color=color, title=title, labels=mid, mae=mae, mape=mapbe)
