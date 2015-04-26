import os
import sys
import glob
import math
import argparse
import itertools
import collections
sys.path.append('/Users/loriab/linux/qcdb')
sys.path.append('/Users/loriab/linux/qcdb/databases')
import qcdb
from qcdb.psivarrosetta import useme2psivar, optclue2psivar
from qcdb.modelchems import Method, BasisSet, methods, bases
from qcdb.exceptions import *
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', 200)


# <<< setup >>>

parser = argparse.ArgumentParser(description='Expand quantum chemical terms into reaction quantities.')
parser.add_argument('-d', '--dbmodule', help='force choice of database module, defaults to project value')
parser.add_argument('-p', '--project', help='select pre-configured project')
parser.add_argument('-v', '--verbose', type=int, default=1, help='amount of printing')
parser.add_argument('-o', '--outdir', default='.', help='directory to write output files')
args = parser.parse_args()

project = args.project
homewrite = args.outdir
verbose = args.verbose

if verbose > 1:
    print args

# <<< read usemefiles and convert to giant DataFrame >>>

#path = r"""C:\Users\Owner\Documents\f12dilabiousemefiles\usemefiles"""
#verbose = True
#homewrite = '/Users/loriab/'
#homewrite = '/Users/loriab/linux/qcdb/sandbox/bfdb'
#homewrite = '/Users/loriab/linux/qcdb/sandbox/flow'
#homewrite = '/Users/loriab/linux/qcdb/sandbox/refitting_dashd'
##homewrite = '/Users/loriab/linux/qcdb/sandbox/ssapt0_for_slipchenko'
##homewrite = '/Users/loriab/linux/qcdb/sandbox/f12dilabio_workup'

# pass --> pass
if project == 'parenq':
    dbse = 'A24'
    ##path = r"""/Users/loriab/linux/qcdb/data/parenqusemefiles/"""
    path = r"""/Users/loriab/linux/qcdb/data/newparenq/"""

# fail --> pass
elif project == 'dft':
    dbse = 'S22'
    path = r"""/Users/loriab/linux/qcdb/data/dftusemefiles/"""

# pass --> pass
elif project == 'f12dilabio':
    dbse = 'A24'
    path = r"""/Users/loriab/linux/qcdb/data/f12dilabiousemefiles/"""

# fail --> pass
elif project == 'dhdft':
    dbse = 'S22'
    path = r"""/Users/loriab/linux/qcdb/data/dhdftusemefiles/"""

#dbse = 'NBC10'
#project == 'dhdft2':
#path = r"""/Users/loriab/linux/qcdb/data/dhdftusemefiles/DSD/"""

# pass --> pass
elif project == 'dilabio':
    dbse = 'A24'
    path = r"""/Users/loriab/linux/qcdb/data/dilabiousemefiles/"""

# pass --> pass
elif project == 'pt2':
    dbse = 'S22'
    path = r"""/Users/loriab/linux/qcdb/data/pt2usemefiles/"""

# pass --> pass
elif project == 'saptone':
    dbse = 'S22'
    path = r"""/Users/loriab/linux/qcdb/data/pt2usemefiles/"""

# not tested
#dbse = 'BBI'
#dbse = 'SSI'
#project = 'merz3'
#path = r"""/Users/loriab/linux/qcdb/data/merz3usemefiles/"""

# not tested
#dbse = 'S22'
#project = 'aep'
#path = r"""/Users/loriab/linux/qcdb/data/pt2usemefiles/"""

elif project == 'sflow':
    dbse = 'S22'
    path = r"""/Users/loriab/linux/qcdb/data/flowusemefiles/"""

elif project == 'dfit':
    dbse = 'HBC6'
    #dbse = 'NBC10ext'
    #dbse = 'ACONF'
    path = r"""/Users/loriab/linux/Refitting_DFT_D/Databases/usemefiles/"""

elif project == 'nbcref':
    dbse = 'NBC10ext'
    path = r"""/Users/loriab/linux/Refitting_DFT_D/Databases/usemefiles/"""

elif project == 'curveref':
    dbse = 'S22by7'
    path = r"""/Users/loriab/linux/Refitting_DFT_D/Databases/usemefiles/"""

else:
    raise ValidationError("""Project %s not defined.""" % (project))

dbse = dbse if args.dbmodule is None else args.dbmodule
dbobj = qcdb.Database(dbse, loadfrompickle=True)
dbse = dbobj.dbse
print '<<<', dbse, '>>>'
modes = []
stoich = {}
for orxn in dbobj.hrxn.itervalues():
    tdict = {}
    for mode, rxnm in orxn.rxnm.iteritems():
        modes.append(mode)
        rgtindx = 0
        for wt in rxnm.itervalues():
            tdict[(mode, 'Rgt' + str(rgtindx))] = wt
            rgtindx += 1
    stoich[orxn.dbrxn] = tdict
dfstoich = pd.DataFrame.from_dict(stoich, orient='index')
modes = set(modes)
names = ['rxn']
names.extend(['Rgt' + str(i) for i in range(12)])
maxrgt = 0
for mode in modes:
    maxrgt = max(maxrgt, dfstoich[mode].shape[1])
print names[:maxrgt+1]
dfstoich.index.names = ['rxn']
h2kc = qcdb.psi_hartree2kcalmol

rawdata = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(dict)))
usemeglob = glob.glob('%s/%s*useme*' % (path, dbse))
if len(usemeglob) == 0:
    raise ValidationError("""No %s usemefiles in %s.""" % (dbse, path))
for useme in usemeglob:
    spl = os.path.basename(useme).split('.')
    basis = spl[0].split('-')[-1]
    piece = '.'.join(spl[1:])
    optns = spl[0].split('-')[2:-1]
    cpmode = 'default'
    optns2 = []
    for opt in sorted(optns):
        if opt == 'CP':
            cpmode = opt
        elif opt == 'unCP':
            cpmode = 'default'
        else:
            try:
                if useme2psivar[piece] in optclue2psivar[opt]:
                    optns2.append(opt)
            except KeyError as e:
                raise ValidationError('Error: option %s needs adding to to optclue2psivar is psivarrosetta.py: %s' % (opt, e))
    optns = '-'.join(optns2)
    #print useme, basis, piece, optns, useme2psivar[piece]
    try:
        print 'Reading useme: %6s %50s %15s %5s' % (basis, useme2psivar[piece], optns, cpmode)
    except KeyError as e:
        raise ValidationError('Error: useme %s needs adding to useme2psivar in psivarrosetta.py' % (e))
        
    if piece.endswith('usemedash'):
        tmp = pd.read_csv('%s' % (useme), index_col=0, sep='\s+', comment='#', na_values='None', names=names[:maxrgt+1])
        rawdata[basis][useme2psivar[piece]][optns]['default'] = tmp.dropna(how='all')
        rawdata[basis][useme2psivar[piece]][optns]['CP'] = tmp.dropna(how='all')
    elif piece.endswith('usemesapt') or piece.endswith('usemedftsapt') or piece.endswith('usemempsapt'):
        # moved labels to top, removed comment marker for labels line, col relabeled to mp2cDisp20 for mpsapt
        tmp = pd.read_csv('%s' % (useme), index_col=0, sep='\s+', comment='#', na_values='None')
        sapt_cols = tmp.columns.tolist()
        tmp = pd.DataFrame(tmp.stack(), columns=['dimer'])
        tmp['dimer'] *= 0.001  # useme in mHartree
        tmp['monoA'] = 0.0
        tmp['monoB'] = 0.0
        tmp = tmp.reorder_levels([1, 0])
        tmp.index.names = ['psivar', 'rxn']
        for pv in sapt_cols:
            try:
                useme2psivar[pv]
            except KeyError as e:
                pass  # bypass extra columns
            else:
                try:
                    tmp2 = tmp.xs(pv, level='psivar')
                except KeyError as e:
                    pass  # bypass empty columns
                else:
                    rawdata[basis][useme2psivar[pv]][optns]['unCP'] = tmp2
                    rawdata[basis][useme2psivar[pv]][optns]['CP'] = tmp2
    else:
        tmp = pd.read_csv('%s' % (useme), index_col=0, sep='\s+', comment='#', na_values='None', names=names[:maxrgt+1])
        rawdata[basis][useme2psivar[piece]][optns][cpmode] = tmp.dropna(how='all')
    #print tmp.head(4)

baszip = {}
for baskey, basval in sorted(rawdata.iteritems()):
    pvzip = {}
    for pvkey, pvval in sorted(basval.iteritems()):
        metazip = {}
        if verbose > 0:
            print """reading  %s %s""" % (pvkey, '.' * (40 - len(pvkey))),
        for metakey, metaval, in pvval.iteritems():
            metazip[metakey] = pd.concat(metaval, axis=1)  # merge CP/unCP
        pvzip[pvkey] = pd.concat(metazip)
        if verbose > 0:
            print """SUCCESS w/ %s %s""" % (baskey, ', '.join(metazip.keys()))
    baszip[baskey] = pd.concat(pvzip)
df = pd.concat(baszip)
df.index.names = ['bstrt', 'psivar', 'meta', 'rxn']


# <<< define utility functions >>>

def rxnm_contract_expand(rgts):
    """Applies the stoichiometry array to the pd.DataFrame *rgts*, then sums 
    across the reagents in each ACTV mode separately thus computing a reaction 
    quantity. The reaction quantity is copied to every reagent within each ACTV 
    mode for ease of further computation and returned.

    """
    stoich_scaled_rgts = dfstoich.mul(rgts)
    micols = stoich_scaled_rgts.columns.values
    modes, rgts = zip(*micols)
    modezip = {}
    for mode in set(modes):
        rgts = [rg for md, rg in micols if md == mode]
        contracted_rxn = stoich_scaled_rgts[mode].sum(axis=1)
        rgtzip = {rg: contracted_rxn for rg in rgts}
        modezip[mode] = pd.concat(rgtzip, axis=1)
    expanded_rgts = pd.concat(modezip, axis=1)
    return expanded_rgts

def reactionate(cpmode, rgts):
    hjkl = pd.DataFrame(dfstoich[cpmode] * rgts[cpmode])
    return hjkl.sum(axis=1)

def ie2(rgts):
    cpmode = 'CP'
    hjkl = pd.DataFrame(dfstoich[cpmode] * rgts[cpmode])
    return hjkl.sum(axis=1)

def ie_uncp2(rgts):
    hjkl = pd.DataFrame(dfstoich['default'] * rgts['default'])
    #hjkl = pd.DataFrame(dfstoich['default'] * rgts['unCP'])
    return hjkl.sum(axis=1)

def default(rgts):
    cpmode = 'default'
    hjkl = pd.DataFrame(dfstoich[cpmode] * rgts[cpmode])
    return hjkl.sum(axis=1)

#def ie(rgts):
#    cpmode = 'CP'
#    return rgts[cpmode]['dimer'] - rgts[cpmode]['monoA'] - rgts[cpmode]['monoB']
#
#
#def ie_uncp(rgts):
#    cpmode = 'unCP'
#    return rgts[cpmode]['dimer'] - rgts[cpmode]['monoA'] - rgts[cpmode]['monoB']
#
#
#def ie_ave(rgts):
#    return 0.5 * (rgts['CP']['dimer'] - rgts['CP']['monoA'] - rgts['CP']['monoB'] +
#                  rgts['unCP']['dimer'] - rgts['unCP']['monoA'] - rgts['unCP']['monoB'])


def categories(df, lvl):
    return sorted(set([tup[lvl] for tup in df.index.values]))


def append_result_of_func_with_funcargs_to_master_DataFrame_atlevel_with_label(master, label, atlevel, func, funcargs):
    odr = {'psivar': [1, 3, 0, 2],
           'bstrt': [3, 1, 0, 2]}
    data_rich_args = [master.xs(pv, level=atlevel) if isinstance(pv, basestring) else pv for pv in funcargs]
    multiopt = []
    for item in data_rich_args:
        try:
            ####print '  C:', sorted(set([tup[1] for tup in item.index.values]))
            multiopt.append(sorted(set([tup[1] for tup in item.index.values])))
            #multiopt.append(set([tup[1] for tup in item.index.values]))
            #multiopt.append(set(['-'.join(sorted(set(tup[1].split('-')))) for tup in item.index.values]))
            #multiopt.append(set(['-'.join(sorted(set(bn.split('-')))) for bn in sorted(set([tup[1] for tup in item.index.values]))]))
            #if label == 'MP2-F12 CORRELATION ENERGY':
            #if label == 'adtz':
            #    bnbn = sorted(set([tup[1] for tup in item.index.values]))
            #    print '    s1:', bnbn
            #    bnbn2 = set(['-'.join(sorted(set(bn.split('-')))) for bn in bnbn])
            #    print '    s2:', bnbn2
                #bnbn3 = 
                #print '    s3:', bnbn3
            #    print '   tup:', ['-'.join(sorted(set(tup[1].split('-')))) for tup in item.index.values]
                #print '   set:', set([tup[1] for tup in item.index.values])
        except AttributeError:
            pass
    #if label == 'MP2-F12 CORRELATION ENERGY':
    #if label == 'adtz':
    ####print '    ... opt comb:', multiopt
    #    print '    ... options combinations:', list(itertools.product(*multiopt))
    optzip = {}
    for multiopt in itertools.product(*multiopt):
        imultiopt = iter(multiopt)
        optlabel = '-'.join(set(sorted([opt for opt in multiopt if opt != ''])))

#        mm = optlabel.split('-')
#        if len(mm) > len(set(mm)):
#            #print '        die nan'
#            continue
#        print '      L:', optlabel, 'M:', multiopt
        data_rich_args_2 = [arg.xs(imultiopt.next(), level='meta') if isinstance(arg, pd.DataFrame) else arg for arg in data_rich_args]
        #if label.startswith('SCS(MI)-MP2-F12') or label == 'adtz':
        #    print label, optlabel
        #    print data_rich_args_2
        
        candidate = func(data_rich_args_2)
        candidate = candidate.dropna(axis=0, how='all')
        if len(candidate) > 0:
            optzip[optlabel] = candidate
            ####print '      L:', optlabel, 'M:', multiopt
        #optzip[optlabel] = func(data_rich_args_2)
    #if label == 'SCS(MI)-MP2-F12 CORRELATION ENERGY':
    #    print 'optzip', optzip
    temp = pd.concat(optzip)
    temp[atlevel] = label
    temp.set_index(atlevel, append=True, inplace=True)
    #print atlevel, 'QQpre', temp.index.values[0]
    temp = temp.reorder_levels(odr[atlevel])
    #print atlevel, 'WWpst', temp.index.values[0]
    return master.append(temp, verify_integrity=True)
    #return master.combine_first(temp)


# <<< define simple functions relating psi variables; follow args use in omega >>>

def xtpl_power(args):
    power, zHI, eHI, eLO = args
    return (eHI * zHI ** power - eLO * (zHI - 1) ** power) / (zHI ** power - (zHI - 1) ** power)


def dispersion_weighting(args):
    omega, hb_mtd, dd_mtd = args
    return omega * hb_mtd + (1.0 - omega) * dd_mtd


def omega(args):
    alpha, beta, ratio = args
    return 0.5 * (1.0 + np.tanh(alpha + beta * ratio))


def difference(args):
    minuend, subtrahend = args
    return minuend - subtrahend


def product(args):
    multiplicand, multiplier = args
    return multiplicand * multiplier


def spin_component_scaling(args):
    os_scale, ss_scale, tot_corl, ss_corl = args
    return os_scale * (tot_corl - ss_corl) + ss_scale * ss_corl
    
    
def correction(args):
    base, plus, minus = args
    return base + plus - minus
    

# <<< append to main DataFrame basic psivar equalities not explicit to useme structure >>>

lvl = 'psivar'
pv0 = collections.OrderedDict()
pv0['MP2-F12 TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv0['CCSD-F12A TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD-F12A CORRELATION ENERGY']}
pv0['CCSD-F12B TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD-F12B CORRELATION ENERGY']}
pv0['CCSD-F12C TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD-F12C CORRELATION ENERGY']}
#pv0['CCSD(T**)-F12A TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T**)-F12A CORRELATION ENERGY']}
#pv0['CCSD(T**)-F12B TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T**)-F12B CORRELATION ENERGY']}
#pv0['CCSD(T**)-F12C TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T**)-F12C CORRELATION ENERGY']}
pv0['CCSD(T)-F12A CORRELATION ENERGY'] = {'func': sum, 'args': ['CCSD-F12A CORRELATION ENERGY', '(T)-F12AB CORRECTION ENERGY']}
pv0['CCSD(T)-F12B CORRELATION ENERGY'] = {'func': sum, 'args': ['CCSD-F12B CORRELATION ENERGY', '(T)-F12AB CORRECTION ENERGY']}
pv0['CCSD(T)-F12C CORRELATION ENERGY'] = {'func': sum, 'args': ['CCSD-F12C CORRELATION ENERGY', '(T)-F12C CORRECTION ENERGY']}
pv0['CCSD(T)-F12A TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T)-F12A CORRELATION ENERGY']}
pv0['CCSD(T)-F12B TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T)-F12B CORRELATION ENERGY']}
pv0['CCSD(T)-F12C TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T)-F12C CORRELATION ENERGY']}
pv0['B3LYP TOTAL ENERGY'] = {'func': sum, 'args': ['B3LYP FUNCTIONAL TOTAL ENERGY']}
pv0['B970 TOTAL ENERGY'] = {'func': sum, 'args': ['B970 FUNCTIONAL TOTAL ENERGY']}
pv0['B97 TOTAL ENERGY'] = {'func': sum, 'args': ['B97 FUNCTIONAL TOTAL ENERGY']}
pv0['BP86 TOTAL ENERGY'] = {'func': sum, 'args': ['BP86 FUNCTIONAL TOTAL ENERGY']}
pv0['BLYP TOTAL ENERGY'] = {'func': sum, 'args': ['BLYP FUNCTIONAL TOTAL ENERGY']}
pv0['M05-2X TOTAL ENERGY'] = {'func': sum, 'args': ['M05-2X FUNCTIONAL TOTAL ENERGY']}
pv0['M06-2X TOTAL ENERGY'] = {'func': sum, 'args': ['M06-2X FUNCTIONAL TOTAL ENERGY']}
pv0['PBE0 TOTAL ENERGY'] = {'func': sum, 'args': ['PBE0 FUNCTIONAL TOTAL ENERGY']}
pv0['PBE TOTAL ENERGY'] = {'func': sum, 'args': ['PBE FUNCTIONAL TOTAL ENERGY']}
pv0['WPBE TOTAL ENERGY'] = {'func': sum, 'args': ['WPBE FUNCTIONAL TOTAL ENERGY']}
pv0['HF TOTAL ENERGY'] = {'func': sum, 'args': ['SCF TOTAL ENERGY']}
pv0['MP3 CORRELATION ENERGY'] = {'func': difference, 'args': ['MP3 FNO CORRELATION ENERGY', 'FNO CORRECTION ENERGY']}  # TODO not checked
pv0['CCSD CORRELATION ENERGY'] = {'func': difference, 'args': ['CCSD FNO CORRELATION ENERGY', 'FNO CORRECTION ENERGY']}  # TODO not checked
pv0['CCSD(T) CORRELATION ENERGY'] = {'func': difference, 'args': ['CCSD(T) FNO CORRELATION ENERGY', 'FNO CORRECTION ENERGY']}
pv0['CCSD(T) TOTAL ENERGY'] = {'func': sum, 'args': ['HF TOTAL ENERGY', 'CCSD(T) CORRELATION ENERGY']}
pv0['CCSDT TOTAL ENERGY'] = {'func': sum, 'args': ['HF TOTAL ENERGY', 'CCSDT CORRELATION ENERGY']}
pv0['CCSDT(Q) TOTAL ENERGY'] = {'func': sum, 'args': ['HF TOTAL ENERGY', 'CCSDT(Q) CORRELATION ENERGY']}
pv0['VV10 TOTAL ENERGY'] = {'func': sum, 'args': ['VV10 FUNCTIONAL TOTAL ENERGY']}
pv0['LC-VV10 TOTAL ENERGY'] = {'func': sum, 'args': ['LC-VV10 FUNCTIONAL TOTAL ENERGY']}
pv0['M08-HX TOTAL ENERGY'] = {'func': sum, 'args': ['M08-HX FUNCTIONAL TOTAL ENERGY']}
pv0['M08-SO TOTAL ENERGY'] = {'func': sum, 'args': ['M08-SO FUNCTIONAL TOTAL ENERGY']}
pv0['M11 TOTAL ENERGY'] = {'func': sum, 'args': ['M11 FUNCTIONAL TOTAL ENERGY']}
pv0['M11L TOTAL ENERGY'] = {'func': sum, 'args': ['M11L FUNCTIONAL TOTAL ENERGY']}
pv0['MP2 TOTAL ENERGY'] = {'func': sum, 'args': ['HF TOTAL ENERGY', 'MP2 CORRELATION ENERGY']}
pv0['CCSD TOTAL ENERGY'] = {'func': sum, 'args': ['HF TOTAL ENERGY', 'CCSD CORRELATION ENERGY']}
pv0['MP3 TOTAL ENERGY'] = {'func': sum, 'args': ['HF TOTAL ENERGY', 'MP3 CORRELATION ENERGY']}
for pvar, action in pv0.iteritems():
    try:
        if verbose > 0:
            print """building %s %s""" % (pvar, '.' * (50 - len(pvar))),
        df = append_result_of_func_with_funcargs_to_master_DataFrame_atlevel_with_label(
            master=df, label=pvar, atlevel=lvl, func=action['func'], funcargs=action['args'])
        if verbose > 0:
            print """SUCCESS"""
    except KeyError, e:
        if verbose > 0:
            print """FAILED, missing %s""" % (e)

# <<< miscellaneous pre-computing >>>
    
#    <<< SAPT EXCHSCAL >>>

try:
    ex10 = df.xs('SAPT EXCH10 ENERGY', level='psivar')
    ex10ss = df.xs('SAPT EXCH10(S^2) ENERGY', level='psivar')
except KeyError, e:
    print 'Not handled: SAPT EXCHSCAL', e
else:
    ex10ss.loc[:,('CP','monoA')] = np.nan
    ex10ss.loc[:,('CP','monoB')] = np.nan
    ex10ss.loc[:,('unCP','monoA')] = np.nan
    ex10ss.loc[:,('unCP','monoB')] = np.nan
    ratio = ex10 / ex10ss
    ones = ratio.copy()
    ones['dimer'] = 1.0
    logic = ex10 > 1.0e-5
    # return 1.0 if ex10 < 1.0e-5 else ex10 / ex10ss
    exsc = ratio.where(logic, ones)
    exsc.loc[:,('CP','monoA')] = 0.0
    exsc.loc[:,('CP','monoB')] = 0.0
    exsc.loc[:,('unCP','monoA')] = 0.0
    exsc.loc[:,('unCP','monoB')] = 0.0
    exsc['psivar'] = 'SAPT EXCHSCAL'
    exsc.set_index('psivar', append=True, inplace=True)
    exsc = exsc.reorder_levels([0, 3, 1, 2])
    exsc.index.names = ['bstrt', 'psivar', 'meta', 'rxn']
    print 'building SAPT EXCHSCAL'
    df = df.append(exsc, verify_integrity=True)

#     <<< DW-MP2 & DW-MP2-F12 >>>

if verbose > 0:
    print 'building intermediates DW-MP2 OMEGA, DW-MP2-F12-OMEGA ...',
dwmp2 = collections.defaultdict(dict)  #lambda: collections.defaultdict(dict))
meta_combos = {
    '': ['', ''],
    'dfmp': ['dfhf', 'dfhf-dfmp']}  # TODO comprehensive? right?

for mt, mtl in meta_combos.iteritems():
    try:
        dwmp2['DW-MP2 OMEGA'][mt] = \
            rxnm_contract_expand(df.xs('HF TOTAL ENERGY', level='psivar').xs(mtl[0], level='meta')) / \
            rxnm_contract_expand(df.xs('MP2 TOTAL ENERGY', level='psivar').xs(mtl[1], level='meta'))
    except KeyError, e:
        pass
    try:
        dwmp2['DW-MP2-F12 OMEGA'][mt] = \
            rxnm_contract_expand(df.xs('HF-CABS TOTAL ENERGY', level='psivar').xs(mtl[0], level='meta')) / \
            rxnm_contract_expand(df.xs('MP2-F12 TOTAL ENERGY', level='psivar').xs(mtl[1], level='meta'))
    except KeyError, e:
        pass

if dwmp2:
    pvzip = {}
    for pvkey, pvval in dwmp2.iteritems():
        metazip = {}
        for metakey, metaval, in pvval.iteritems():
            metazip[metakey] = metaval                
        pvzip[pvkey] = pd.concat(metazip)
    df_omega = pd.concat(pvzip)

    df_omega = df_omega.reorder_levels([2, 0, 1, 3])
    df_omega.index.names = ['bstrt', 'psivar', 'meta', 'rxn']
    try:
        df_omega = omega([0.15276, 1.89952, df_omega])
    except KeyError, e:
        pass
    if verbose > 0:
        print 'SUCCESS'
    df = df.append(df_omega, verify_integrity=True)
else:
    if verbose > 0:
        print 'Not Handled'

#     <<< DW-CCSD(T)-F12 >>>

try:
    df.xs('HF-CABS TOTAL ENERGY', level='psivar')
    df.xs('MP2-F12 TOTAL ENERGY', level='psivar')
except KeyError, e:
    print 'Not handled: DW-CCSD(T) OMEGA', e
else:
    # LAB 23 Apr 2015 switching from ie
    #print df.xs('HF-CABS TOTAL ENERGY', level='psivar').to_dict()
    ratio_HFCABS_MP2F12 = rxnm_contract_expand(df.xs('HF-CABS TOTAL ENERGY', level='psivar')) / \
                          rxnm_contract_expand(df.xs('MP2-F12 TOTAL ENERGY', level='psivar'))
    print 'ratio_HFCABS_MP2F12\n', ratio_HFCABS_MP2F12
    # TODO possibly set SA or other rxnm to NaN as inapplicable
    dwcc = {}
    try:
        dwcc['adz'] = omega([-1, 4, ratio_HFCABS_MP2F12.xs('adz', level='bstrt')])
        print 'post dwcc adz'
    except KeyError as e:
        pass
    try:
        dwcc['atz'] = omega([0.4, 0.6, ratio_HFCABS_MP2F12.xs('atz', level='bstrt')])
        print 'post dwcc atz'
    except KeyError as e:
        pass
    df_omega = pd.concat(dwcc)
    print 'DFOM#GA0\n', df_omega
    df_omega['psivar'] = 'DW-CCSD(T)-F12 OMEGA'
    df_omega.set_index('psivar', append=True, inplace=True)
    print 'DFOM#GA1\n', df_omega
    df_omega = df_omega.reorder_levels([0, 3, 1, 2])
    df_omega.index.names = ['bstrt', 'psivar', 'meta', 'rxn']
    print 'DFOM#GA2\n', df_omega
    df = df.append(df_omega, verify_integrity=True)

#     <<< (T*)-F12 & (T**)-F12 >>>
try:
    df.xs('MP2 CORRELATION ENERGY', level='psivar')
    df.xs('MP2-F12 CORRELATION ENERGY', level='psivar')
except KeyError, e:
    print 'Not handled: (T*)-F12 & (T**)-F12 SCALE', e
else:
    ratio_MP2F12_MP2 = df.xs('MP2-F12 CORRELATION ENERGY', level='psivar') / df.xs('MP2 CORRELATION ENERGY', level='psivar')
    tstar = ratio_MP2F12_MP2.copy()
    tstar['psivar'] = '(T*)-F12 SCALE'
    tstar.set_index('psivar', append=True, inplace=True)
    tstar = tstar.reorder_levels([0, 3, 1, 2])
    tstar.index.names = ['bstrt', 'psivar', 'meta', 'rxn']
    df = df.append(tstar, verify_integrity=True)
    print 'post tstart'
        
    tstarstar = ratio_MP2F12_MP2.copy()
    for md, rg in tstarstar.columns.values:
        if rg != 'Rgt0':
            tstarstar.loc[:, (md, rg)] = tstarstar[md]['Rgt0']
    tstarstar['psivar'] = '(T**)-F12 SCALE'
    tstarstar.set_index('psivar', append=True, inplace=True)
    tstarstar = tstarstar.reorder_levels([0, 3, 1, 2])
    tstarstar.index.names = ['bstrt', 'psivar', 'meta', 'rxn']
    print 'post tsrarstar'
    df = df.append(tstarstar, verify_integrity=True)

#     <<< DASH-D >>>
try:
    df.xs('nobas', level='bstrt')
except KeyError, e:
    print 'Not handled: DASH-D', e
else:
    print 'indash'
    df_nobas = pd.concat({tup[0]: df.xs('nobas', level='bstrt') for tup in df.index.values if tup[0] != 'nobas'})
    print df_nobas.head(50)
    df_nobas.index.names = ['bstrt', 'psivar', 'meta', 'rxn']
    print 'postdash'
    df = df.append(df_nobas, verify_integrity=True)


#mlist = ['B3LYP FUNCTIONAL TOTAL ENERGY', 'B3LYP-D2 DISPERSION CORRECTION ENERGY']
#for ml in mlist:
#    mdf = df.xs(ml, level='psivar').xs('adz', level='bstrt')
#    #mdf = df.xs(ml, level='psivar').xs('adz', level='bstrt').xs('S22-2', level='rxn')
#    print ml, '\n'#, categories(mdf, 0)
#    #print h2kc * mdf
#    print mdf.head(5)

# <<< append to main DataFrame computable method quantities >>>

lvl = 'psivar'
pv1 = collections.OrderedDict()
pv1['B3LYP-D2 TOTAL ENERGY'] = {'func': sum, 'args': ['B3LYP FUNCTIONAL TOTAL ENERGY', 'B3LYP-D2 DISPERSION CORRECTION ENERGY']}
pv1['B3LYP-D3 TOTAL ENERGY'] = {'func': sum, 'args': ['B3LYP FUNCTIONAL TOTAL ENERGY', 'B3LYP-D3 DISPERSION CORRECTION ENERGY']}
pv1['B3LYP-D3BJ TOTAL ENERGY'] = {'func': sum, 'args': ['B3LYP FUNCTIONAL TOTAL ENERGY', 'B3LYP-D3BJ DISPERSION CORRECTION ENERGY']}
pv1['B3LYP-XDM TOTAL ENERGY'] = {'func': sum, 'args': ['B3LYP FUNCTIONAL TOTAL ENERGY', 'B3LYP-XDM DISPERSION CORRECTION ENERGY']}
pv1['B2PLYP-D2 TOTAL ENERGY'] = {'func': sum, 'args': ['B2PLYP TOTAL ENERGY', 'B2PLYP-D2 DISPERSION CORRECTION ENERGY']}
pv1['B2PLYP-D3 TOTAL ENERGY'] = {'func': sum, 'args': ['B2PLYP TOTAL ENERGY', 'B2PLYP-D3 DISPERSION CORRECTION ENERGY']}
pv1['B2PLYP-D3BJ TOTAL ENERGY'] = {'func': sum, 'args': ['B2PLYP TOTAL ENERGY', 'B2PLYP-D3BJ DISPERSION CORRECTION ENERGY']}
pv1['DSD-PBEP86-D2 TOTAL ENERGY'] = {'func': sum, 'args': ['DSD-PBEP86 TOTAL ENERGY', 'DSD-PBEP86-D2 DISPERSION CORRECTION ENERGY']}
pv1['DSD-PBEP86-D3 TOTAL ENERGY'] = {'func': sum, 'args': ['DSD-PBEP86 TOTAL ENERGY', 'DSD-PBEP86-D3 DISPERSION CORRECTION ENERGY']}
pv1['DSD-PBEP86-D3BJ TOTAL ENERGY'] = {'func': sum, 'args': ['DSD-PBEP86 TOTAL ENERGY', 'DSD-PBEP86-D3BJ DISPERSION CORRECTION ENERGY']}
pv1['B970-D2 TOTAL ENERGY'] = {'func': sum, 'args': ['B970 FUNCTIONAL TOTAL ENERGY', 'B970-D2 DISPERSION CORRECTION ENERGY']}
pv1['B97-D2 TOTAL ENERGY'] = {'func': sum, 'args': ['B97 FUNCTIONAL TOTAL ENERGY', 'B97-D2 DISPERSION CORRECTION ENERGY']}
pv1['B97-D3 TOTAL ENERGY'] = {'func': sum, 'args': ['B97 FUNCTIONAL TOTAL ENERGY', 'B97-D3 DISPERSION CORRECTION ENERGY']}
pv1['B97-D3BJ TOTAL ENERGY'] = {'func': sum, 'args': ['B97 FUNCTIONAL TOTAL ENERGY', 'B97-D3BJ DISPERSION CORRECTION ENERGY']}
pv1['BP86-D2 TOTAL ENERGY'] = {'func': sum, 'args': ['BP86 FUNCTIONAL TOTAL ENERGY', 'BP86-D2 DISPERSION CORRECTION ENERGY']}
pv1['BP86-D3 TOTAL ENERGY'] = {'func': sum, 'args': ['BP86 FUNCTIONAL TOTAL ENERGY', 'BP86-D3 DISPERSION CORRECTION ENERGY']}
pv1['BP86-D3BJ TOTAL ENERGY'] = {'func': sum, 'args': ['BP86 FUNCTIONAL TOTAL ENERGY', 'BP86-D3BJ DISPERSION CORRECTION ENERGY']}
pv1['BLYP-D2 TOTAL ENERGY'] = {'func': sum, 'args': ['BLYP FUNCTIONAL TOTAL ENERGY', 'BLYP-D2 DISPERSION CORRECTION ENERGY']}
pv1['BLYP-D3 TOTAL ENERGY'] = {'func': sum, 'args': ['BLYP FUNCTIONAL TOTAL ENERGY', 'BLYP-D3 DISPERSION CORRECTION ENERGY']}
pv1['BLYP-D3BJ TOTAL ENERGY'] = {'func': sum, 'args': ['BLYP FUNCTIONAL TOTAL ENERGY', 'BLYP-D3BJ DISPERSION CORRECTION ENERGY']}
pv1['PBE-D2 TOTAL ENERGY'] = {'func': sum, 'args': ['PBE FUNCTIONAL TOTAL ENERGY', 'PBE-D2 DISPERSION CORRECTION ENERGY']}
pv1['PBE-D3 TOTAL ENERGY'] = {'func': sum, 'args': ['PBE FUNCTIONAL TOTAL ENERGY', 'PBE-D3 DISPERSION CORRECTION ENERGY']}
pv1['PBE-D3BJ TOTAL ENERGY'] = {'func': sum, 'args': ['PBE FUNCTIONAL TOTAL ENERGY', 'PBE-D3BJ DISPERSION CORRECTION ENERGY']}
pv1['PBE0-D2 TOTAL ENERGY'] = {'func': sum, 'args': ['PBE0 FUNCTIONAL TOTAL ENERGY', 'PBE0-D2 DISPERSION CORRECTION ENERGY']}
pv1['PBE0-D3 TOTAL ENERGY'] = {'func': sum, 'args': ['PBE0 FUNCTIONAL TOTAL ENERGY', 'PBE0-D3 DISPERSION CORRECTION ENERGY']}
pv1['PBE0-D3BJ TOTAL ENERGY'] = {'func': sum, 'args': ['PBE0 FUNCTIONAL TOTAL ENERGY', 'PBE0-D3BJ DISPERSION CORRECTION ENERGY']}
pv1['M05-2X-D3 TOTAL ENERGY'] = {'func': sum, 'args': ['M05-2X FUNCTIONAL TOTAL ENERGY', 'M05-2X-D3 DISPERSION CORRECTION ENERGY']}
pv1['M06-2X-D3 TOTAL ENERGY'] = {'func': sum, 'args': ['M06-2X FUNCTIONAL TOTAL ENERGY', 'M06-2X-D3 DISPERSION CORRECTION ENERGY']}
pv1['WB97X-D TOTAL ENERGY'] = {'func': sum, 'args': ['WB97X FUNCTIONAL TOTAL ENERGY', 'WB97X-D DISPERSION CORRECTION ENERGY']}
pv1['DLDF+D TOTAL ENERGY'] = {'func': sum, 'args': ['DLDF FUNCTIONAL TOTAL ENERGY', 'DLDF+D DISPERSION CORRECTION ENERGY']}
pv1['CCSD(T) CC CORRECTION ENERGY'] = {'func': difference, 'args': ['CCSD(T) CORRELATION ENERGY', 'MP2 CORRELATION ENERGY']}
pv1['MP3 CC CORRECTION ENERGY'] = {'func': difference, 'args': ['MP3 CORRELATION ENERGY', 'MP2 CORRELATION ENERGY']}
pv1['MP2.5 CORRELATION ENERGY'] = {'func': dispersion_weighting, 'args': [0.5, 'MP3 CORRELATION ENERGY', 'MP2 CORRELATION ENERGY']}
pv1['MP2.5 TOTAL ENERGY'] = {'func': sum, 'args': ['HF TOTAL ENERGY', 'MP2.5 CORRELATION ENERGY']}
pv1['MP2.5 CC CORRECTION ENERGY'] = {'func': difference, 'args': ['MP2.5 CORRELATION ENERGY', 'MP2 CORRELATION ENERGY']}
pv1['CCSD CC CORRECTION ENERGY'] = {'func': difference, 'args': ['CCSD CORRELATION ENERGY', 'MP2 CORRELATION ENERGY']}
pv1['SCS-CCSD CC CORRECTION ENERGY'] = {'func': difference, 'args': ['SCS-CCSD CORRELATION ENERGY', 'MP2 CORRELATION ENERGY']}
pv1['SCS(MI)-CCSD CC CORRECTION ENERGY'] = {'func': difference, 'args': ['SCS(MI)-CCSD CORRELATION ENERGY', 'MP2 CORRELATION ENERGY']}
pv1['SCS-MP2 CORRELATION ENERGY'] = {'func': spin_component_scaling, 'args': [1.2, 1.0/3.0, 'MP2 CORRELATION ENERGY', 'MP2 SAME-SPIN CORRELATION ENERGY']}
pv1['SCS-MP2 TOTAL ENERGY'] = {'func': sum, 'args': ['HF TOTAL ENERGY', 'SCS-MP2 CORRELATION ENERGY']}
pv1['SCS(N)-MP2 CORRELATION ENERGY'] = {'func': spin_component_scaling, 'args': [0.0, 1.76, 'MP2 CORRELATION ENERGY', 'MP2 SAME-SPIN CORRELATION ENERGY']}
pv1['SCS(N)-MP2 TOTAL ENERGY'] = {'func': sum, 'args': ['HF TOTAL ENERGY', 'SCS(N)-MP2 CORRELATION ENERGY']}
pv1['DW-MP2 CORRELATION ENERGY'] = {'func': dispersion_weighting, 'args': ['DW-MP2 OMEGA', 'MP2 CORRELATION ENERGY', 'SCS-MP2 CORRELATION ENERGY']}
pv1['DW-MP2 TOTAL ENERGY'] = {'func': sum, 'args': ['HF TOTAL ENERGY', 'DW-MP2 CORRELATION ENERGY']}
pv1['SCS-MP2-F12 CORRELATION ENERGY'] = {'func': spin_component_scaling, 'args': [1.2, 1.0/3.0, 'MP2-F12 CORRELATION ENERGY', 'MP2-F12 SAME-SPIN CORRELATION ENERGY']}
pv1['SCS-MP2-F12 TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'SCS-MP2-F12 CORRELATION ENERGY']}
pv1['SCS(N)-MP2-F12 CORRELATION ENERGY'] = {'func': spin_component_scaling, 'args': [0.0, 1.76, 'MP2-F12 CORRELATION ENERGY', 'MP2-F12 SAME-SPIN CORRELATION ENERGY']}
pv1['SCS(N)-MP2-F12 TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'SCS(N)-MP2-F12 CORRELATION ENERGY']}
pv1['DW-MP2-F12 CORRELATION ENERGY'] = {'func': dispersion_weighting, 'args': ['DW-MP2-F12 OMEGA', 'MP2-F12 CORRELATION ENERGY', 'SCS-MP2-F12 CORRELATION ENERGY']}
pv1['DW-MP2-F12 TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'DW-MP2-F12 CORRELATION ENERGY']}
pv1['SCS-CCSD-F12A CORRELATION ENERGY'] = {'func': spin_component_scaling, 'args': [1.27, 1.13, 'CCSD-F12A CORRELATION ENERGY', 'CCSD-F12A SAME-SPIN CORRELATION ENERGY']}
pv1['SCS-CCSD-F12A TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'SCS-CCSD-F12A CORRELATION ENERGY']}
pv1['SCS(MI)-CCSD-F12A CORRELATION ENERGY'] = {'func': spin_component_scaling, 'args': [1.11, 1.28, 'CCSD-F12A CORRELATION ENERGY', 'CCSD-F12A SAME-SPIN CORRELATION ENERGY']}
pv1['SCS(MI)-CCSD-F12A TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'SCS(MI)-CCSD-F12A CORRELATION ENERGY']}
pv1['SCS-CCSD-F12B CORRELATION ENERGY'] = {'func': spin_component_scaling, 'args': [1.27, 1.13, 'CCSD-F12B CORRELATION ENERGY', 'CCSD-F12B SAME-SPIN CORRELATION ENERGY']}
pv1['SCS-CCSD-F12B TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'SCS-CCSD-F12B CORRELATION ENERGY']}
pv1['SCS(MI)-CCSD-F12B CORRELATION ENERGY'] = {'func': spin_component_scaling, 'args': [1.11, 1.28, 'CCSD-F12B CORRELATION ENERGY', 'CCSD-F12B SAME-SPIN CORRELATION ENERGY']}
pv1['SCS(MI)-CCSD-F12B TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'SCS(MI)-CCSD-F12B CORRELATION ENERGY']}
pv1['SCS-CCSD-F12C CORRELATION ENERGY'] = {'func': spin_component_scaling, 'args': [1.27, 1.13, 'CCSD-F12C CORRELATION ENERGY', 'CCSD-F12C SAME-SPIN CORRELATION ENERGY']}
pv1['SCS-CCSD-F12C TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'SCS-CCSD-F12C CORRELATION ENERGY']}
pv1['SCS(MI)-CCSD-F12C CORRELATION ENERGY'] = {'func': spin_component_scaling, 'args': [1.11, 1.28, 'CCSD-F12C CORRELATION ENERGY', 'CCSD-F12C SAME-SPIN CORRELATION ENERGY']}
pv1['SCS(MI)-CCSD-F12C TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'SCS(MI)-CCSD-F12C CORRELATION ENERGY']}
pv1['CCSD-F12A CC CORRECTION ENERGY'] = {'func': difference, 'args': ['CCSD-F12A CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['CCSD-F12B CC CORRECTION ENERGY'] = {'func': difference, 'args': ['CCSD-F12B CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['CCSD-F12C CC CORRECTION ENERGY'] = {'func': difference, 'args': ['CCSD-F12C CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['SCS-CCSD-F12A CC CORRECTION ENERGY'] = {'func': difference, 'args': ['SCS-CCSD-F12A CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['SCS-CCSD-F12B CC CORRECTION ENERGY'] = {'func': difference, 'args': ['SCS-CCSD-F12B CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['SCS-CCSD-F12C CC CORRECTION ENERGY'] = {'func': difference, 'args': ['SCS-CCSD-F12C CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['SCS(MI)-CCSD-F12A CC CORRECTION ENERGY'] = {'func': difference, 'args': ['SCS(MI)-CCSD-F12A CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['SCS(MI)-CCSD-F12B CC CORRECTION ENERGY'] = {'func': difference, 'args': ['SCS(MI)-CCSD-F12B CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['SCS(MI)-CCSD-F12C CC CORRECTION ENERGY'] = {'func': difference, 'args': ['SCS(MI)-CCSD-F12C CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['(T*)-F12AB CORRECTION ENERGY'] = {'func': product, 'args': ['(T*)-F12 SCALE', '(T)-F12AB CORRECTION ENERGY']}
pv1['(T*)-F12C CORRECTION ENERGY'] = {'func': product, 'args': ['(T*)-F12 SCALE', '(T)-F12C CORRECTION ENERGY']}
pv1['(T**)-F12AB CORRECTION ENERGY'] = {'func': product, 'args': ['(T**)-F12 SCALE', '(T)-F12AB CORRECTION ENERGY']}
pv1['(T**)-F12C CORRECTION ENERGY'] = {'func': product, 'args': ['(T**)-F12 SCALE', '(T)-F12C CORRECTION ENERGY']}
pv1['CCSD(T)-F12A CC CORRECTION ENERGY'] = {'func': difference, 'args': ['CCSD(T)-F12A CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['CCSD(T)-F12B CC CORRECTION ENERGY'] = {'func': difference, 'args': ['CCSD(T)-F12B CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['CCSD(T)-F12C CC CORRECTION ENERGY'] = {'func': difference, 'args': ['CCSD(T)-F12C CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['CCSD(T*)-F12A CORRELATION ENERGY'] = {'func': sum, 'args': ['CCSD-F12A CORRELATION ENERGY', '(T*)-F12AB CORRECTION ENERGY']}
pv1['CCSD(T*)-F12B CORRELATION ENERGY'] = {'func': sum, 'args': ['CCSD-F12B CORRELATION ENERGY', '(T*)-F12AB CORRECTION ENERGY']}
pv1['CCSD(T*)-F12C CORRELATION ENERGY'] = {'func': sum, 'args': ['CCSD-F12C CORRELATION ENERGY', '(T*)-F12C CORRECTION ENERGY']}
pv1['CCSD(T*)-F12A CC CORRECTION ENERGY'] = {'func': difference, 'args': ['CCSD(T*)-F12A CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['CCSD(T*)-F12B CC CORRECTION ENERGY'] = {'func': difference, 'args': ['CCSD(T*)-F12B CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['CCSD(T*)-F12C CC CORRECTION ENERGY'] = {'func': difference, 'args': ['CCSD(T*)-F12C CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['CCSD(T*)-F12A TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T*)-F12A CORRELATION ENERGY']}
pv1['CCSD(T*)-F12B TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T*)-F12B CORRELATION ENERGY']}
pv1['CCSD(T*)-F12C TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T*)-F12C CORRELATION ENERGY']}
pv1['CCSD(T**)-F12A CORRELATION ENERGY'] = {'func': sum, 'args': ['CCSD-F12A CORRELATION ENERGY', '(T**)-F12AB CORRECTION ENERGY']}
pv1['CCSD(T**)-F12B CORRELATION ENERGY'] = {'func': sum, 'args': ['CCSD-F12B CORRELATION ENERGY', '(T**)-F12AB CORRECTION ENERGY']}
pv1['CCSD(T**)-F12C CORRELATION ENERGY'] = {'func': sum, 'args': ['CCSD-F12C CORRELATION ENERGY', '(T**)-F12C CORRECTION ENERGY']}
pv1['CCSD(T**)-F12A CC CORRECTION ENERGY'] = {'func': difference, 'args': ['CCSD(T**)-F12A CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['CCSD(T**)-F12B CC CORRECTION ENERGY'] = {'func': difference, 'args': ['CCSD(T**)-F12B CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['CCSD(T**)-F12C CC CORRECTION ENERGY'] = {'func': difference, 'args': ['CCSD(T**)-F12C CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['CCSD(T**)-F12A TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T**)-F12A CORRELATION ENERGY']}  # duplicate def
pv1['CCSD(T**)-F12B TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T**)-F12B CORRELATION ENERGY']}  # duplicate def
pv1['CCSD(T**)-F12C TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T**)-F12C CORRELATION ENERGY']}  # duplicate def
pv1['DW-CCSD(T)-F12 CORRELATION ENERGY'] = {'func': dispersion_weighting, 'args': ['DW-CCSD(T)-F12 OMEGA', 'CCSD(T)-F12A CORRELATION ENERGY', 'CCSD(T)-F12B CORRELATION ENERGY']}
pv1['DW-CCSD(T)-F12 TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'DW-CCSD(T)-F12 CORRELATION ENERGY']}
pv1['DW-CCSD(T)-F12 CC CORRECTION ENERGY'] = {'func': difference, 'args': ['DW-CCSD(T)-F12 CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['DW-CCSD(T*)-F12 CORRELATION ENERGY'] = {'func': dispersion_weighting, 'args': ['DW-CCSD(T)-F12 OMEGA', 'CCSD(T*)-F12A CORRELATION ENERGY', 'CCSD(T*)-F12B CORRELATION ENERGY']}
pv1['DW-CCSD(T*)-F12 TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'DW-CCSD(T*)-F12 CORRELATION ENERGY']}
pv1['DW-CCSD(T*)-F12 CC CORRECTION ENERGY'] = {'func': difference, 'args': ['DW-CCSD(T*)-F12 CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['DW-CCSD(T**)-F12 CORRELATION ENERGY'] = {'func': dispersion_weighting, 'args': ['DW-CCSD(T)-F12 OMEGA', 'CCSD(T**)-F12A CORRELATION ENERGY', 'CCSD(T**)-F12B CORRELATION ENERGY']}
pv1['DW-CCSD(T**)-F12 TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'DW-CCSD(T**)-F12 CORRELATION ENERGY']}
pv1['DW-CCSD(T**)-F12 CC CORRECTION ENERGY'] = {'func': difference, 'args': ['DW-CCSD(T**)-F12 CORRELATION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['MP2C CC CORRECTION ENERGY'] = {'func': difference, 'args': ['MP2C DISP20 ENERGY', 'SAPT DISP20 ENERGY']}
pv1['MP2C CORRELATION ENERGY'] = {'func': sum, 'args': ['MP2C CC CORRECTION ENERGY', 'MP2 CORRELATION ENERGY']}
pv1['MP2C TOTAL ENERGY'] = {'func': sum, 'args': ['HF TOTAL ENERGY', 'MP2C CORRELATION ENERGY']}
pv1['MP2C-F12 CC CORRECTION ENERGY'] = {'func': sum, 'args': ['MP2C CC CORRECTION ENERGY']}
pv1['MP2C-F12 CORRELATION ENERGY'] = {'func': sum, 'args': ['MP2C CC CORRECTION ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv1['MP2C-F12 TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'MP2C-F12 CORRELATION ENERGY']}
pv1['SAPT EXCHSCAL3'] = {'func': lambda x: x[0] ** 3, 'args': ['SAPT EXCHSCAL']}
pv1['SAPT HF(2) ALPHA=0.0 ENERGY'] = {'func': lambda x: x[0] - (x[1] + x[2] + x[3] + x[4]), 
                            'args': ['SAPT HF TOTAL ENERGY', 'SAPT ELST10,R ENERGY', 'SAPT EXCH10 ENERGY', 
                                     'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY']}
pv1['SAPT HF(2) ENERGY'] = {'func': lambda x: x[1] + (1.0 - x[0]) * x[2],
                                      'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ALPHA=0.0 ENERGY', 'SAPT EXCH-IND20,R ENERGY']}
pv1['SAPT HF(3) ENERGY'] = {'func': lambda x: x[1] - (x[2] + x[0] * x[3]),
                            'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ENERGY', 'SAPT IND30,R ENERGY', 'SAPT EXCH-IND30,R ENERGY']}
pv1['SAPT MP2(2) ENERGY'] = {'func': lambda x: x[1] - (x[2] + x[3] + x[4] + x[0] * (x[5] + x[6] + x[7] + x[8])),
                             'args': ['SAPT EXCHSCAL', 'MP2 CORRELATION ENERGY', 'SAPT ELST12,R ENERGY', 
                                      'SAPT IND22 ENERGY', 'SAPT DISP20 ENERGY', 'SAPT EXCH11(S^2) ENERGY', 
                                      'SAPT EXCH12(S^2) ENERGY', 'SAPT EXCH-IND22 ENERGY', 'SAPT EXCH-DISP20 ENERGY']}
pv1['SAPT MP2(3) ENERGY'] = {'func': lambda x: x[1] - (x[2] + x[0] * x[3]),
                             'args': ['SAPT EXCHSCAL', 'SAPT MP2(2) ENERGY', 'SAPT IND-DISP30 ENERGY', 'SAPT EXCH-IND-DISP30 ENERGY']}
pv1['SAPT MP4 DISP'] = {'func': lambda x: x[0] * x[1] + x[2] + x[3] + x[4] + x[5],
                        'args': ['SAPT EXCHSCAL', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP20 ENERGY',
                                 'SAPT DISP21 ENERGY', 'SAPT DISP22(SDQ) ENERGY', 'SAPT EST.DISP22(T) ENERGY']}
pv1['SAPT CCD DISP'] = {'func': lambda x: x[0] * x[1] + x[2] + x[3] + x[4],
                        'args': ['SAPT EXCHSCAL', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP2(CCD) ENERGY', 
                                 'SAPT DISP22(S)(CCD) ENERGY', 'SAPT EST.DISP22(T)(CCD) ENERGY']}
pv1['SAPT0 ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY']}
pv1['SAPT0 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT EXCH10 ENERGY']}
pv1['SAPT0 INDC ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3], 
                            'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ENERGY', 'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY']}
pv1['SAPT0 DISP ENERGY'] = {'func': lambda x: x[0] * x[1] + x[2], 
                            'args': ['SAPT EXCHSCAL', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP20 ENERGY']}
pv1['SAPT0 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT0 ELST ENERGY', 'SAPT0 EXCH ENERGY', 'SAPT0 INDC ENERGY', 'SAPT0 DISP ENERGY']}
pv1['SSAPT0 ELST ENERGY'] = {'func': sum, 'args': ['SAPT0 ELST ENERGY']}
pv1['SSAPT0 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT0 EXCH ENERGY']}
pv1['SSAPT0 INDC ENERGY'] = {'func': lambda x: x[1] + (x[0] - 1.0) * x[2],
                            'args': ['SAPT EXCHSCAL3', 'SAPT0 INDC ENERGY', 'SAPT EXCH-IND20,R ENERGY']}
pv1['SSAPT0 DISP ENERGY'] = {'func': lambda x: x[0] * x[1] + x[2], 
                            'args': ['SAPT EXCHSCAL3', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP20 ENERGY']}
pv1['SSAPT0 TOTAL ENERGY'] = {'func': sum, 'args': ['SSAPT0 ELST ENERGY', 'SSAPT0 EXCH ENERGY', 'SSAPT0 INDC ENERGY', 'SSAPT0 DISP ENERGY']}
pv1['SCS-SAPT0 ELST ENERGY'] = {'func': sum, 'args': ['SAPT0 ELST ENERGY']}
pv1['SCS-SAPT0 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT0 EXCH ENERGY']}
pv1['SCS-SAPT0 INDC ENERGY'] = {'func': sum, 'args': ['SAPT0 INDC ENERGY']}
pv1['SCS-SAPT0 DISP ENERGY'] = {'func': lambda x: x[0] * (x[1] + x[2]) + x[3] * (x[4] + x[5]), 
                              'args': [0.66, 'SAPT EXCH-DISP20(SS) ENERGY', 'SAPT DISP20(SS) ENERGY',
                                       1.2, 'SAPT EXCH-DISP20(OS) ENERGY', 'SAPT DISP20(OS) ENERGY']}  # note no xs for SCS disp
pv1['SCS-SAPT0 TOTAL ENERGY'] = {'func': sum, 'args': ['SCS-SAPT0 ELST ENERGY', 'SCS-SAPT0 EXCH ENERGY', 'SCS-SAPT0 INDC ENERGY', 'SCS-SAPT0 DISP ENERGY']}
pv1['SAPT2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY', 'SAPT ELST12,R ENERGY']}
pv1['SAPT2 EXCH ENERGY'] = {'func': lambda x: x[1] + x[0] * (x[2] + x[3]), 
                            'args': ['SAPT EXCHSCAL', 'SAPT EXCH10 ENERGY', 'SAPT EXCH11(S^2) ENERGY', 'SAPT EXCH12(S^2) ENERGY']}
pv1['SAPT2 INDC ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5], 
                            'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ENERGY', 'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY',
                                     'SAPT IND22 ENERGY', 'SAPT EXCH-IND22 ENERGY']}
pv1['SAPT2 DISP ENERGY'] = {'func': lambda x: x[0] * x[1] + x[2], 
                            'args': ['SAPT EXCHSCAL', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP20 ENERGY']}
pv1['SAPT2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2 ELST ENERGY', 'SAPT2 EXCH ENERGY', 'SAPT2 INDC ENERGY', 'SAPT2 DISP ENERGY']}
pv1['SAPT2+ ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY', 'SAPT ELST12,R ENERGY']}
pv1['SAPT2+ EXCH ENERGY'] = {'func': lambda x: x[1] + x[0] * (x[2] + x[3]), 
                             'args': ['SAPT EXCHSCAL', 'SAPT EXCH10 ENERGY', 'SAPT EXCH11(S^2) ENERGY', 'SAPT EXCH12(S^2) ENERGY']}
pv1['SAPT2+ INDC ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5], 
                             'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ENERGY', 'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY',
                                      'SAPT IND22 ENERGY', 'SAPT EXCH-IND22 ENERGY']}
pv1['SAPT2+ DISP ENERGY'] = {'func': sum, 'args': ['SAPT MP4 DISP']}
pv1['SAPT2+ TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+ ELST ENERGY', 'SAPT2+ EXCH ENERGY', 'SAPT2+ INDC ENERGY', 'SAPT2+ DISP ENERGY']}
pv1['SAPT2+(CCD) ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+ ELST ENERGY']}
pv1['SAPT2+(CCD) EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+ EXCH ENERGY']}
pv1['SAPT2+(CCD) INDC ENERGY'] = {'func': sum, 'args': ['SAPT2+ INDC ENERGY']}
pv1['SAPT2+(CCD) DISP ENERGY'] = {'func': sum, 'args': ['SAPT CCD DISP']}
pv1['SAPT2+(CCD) TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(CCD) ELST ENERGY', 'SAPT2+(CCD) EXCH ENERGY', 'SAPT2+(CCD) INDC ENERGY', 'SAPT2+(CCD) DISP ENERGY']}
pv1['SAPT2+DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+ ELST ENERGY']}
pv1['SAPT2+DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+ EXCH ENERGY']}
pv1['SAPT2+DMP2 INDC ENERGY'] = {'func': sum, 'args': ['SAPT2+ INDC ENERGY', 'SAPT MP2(2) ENERGY']}
pv1['SAPT2+DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+ DISP ENERGY']}
pv1['SAPT2+DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+DMP2 ELST ENERGY', 'SAPT2+DMP2 EXCH ENERGY', 'SAPT2+DMP2 INDC ENERGY', 'SAPT2+DMP2 DISP ENERGY']}
pv1['SAPT2+(CCD)DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+ ELST ENERGY']}
pv1['SAPT2+(CCD)DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+ EXCH ENERGY']}
pv1['SAPT2+(CCD)DMP2 INDC ENERGY'] = {'func': sum, 'args': ['SAPT2+DMP2 INDC ENERGY']}
pv1['SAPT2+(CCD)DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+(CCD) DISP ENERGY']}
pv1['SAPT2+(CCD)DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(CCD)DMP2 ELST ENERGY', 'SAPT2+(CCD)DMP2 EXCH ENERGY', 'SAPT2+(CCD)DMP2 INDC ENERGY', 'SAPT2+(CCD)DMP2 DISP ENERGY']}
pv1['SAPT2+(3) ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY', 'SAPT ELST12,R ENERGY', 'SAPT ELST13,R ENERGY']}
pv1['SAPT2+(3) EXCH ENERGY'] = {'func': lambda x: x[1] + x[0] * (x[2] + x[3]), 
                                'args': ['SAPT EXCHSCAL', 'SAPT EXCH10 ENERGY', 'SAPT EXCH11(S^2) ENERGY', 'SAPT EXCH12(S^2) ENERGY']}
pv1['SAPT2+(3) INDC ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5], 
                                'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ENERGY', 'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY',
                                         'SAPT IND22 ENERGY', 'SAPT EXCH-IND22 ENERGY']}
pv1['SAPT2+(3) DISP ENERGY'] = {'func': sum, 'args': ['SAPT MP4 DISP', 'SAPT DISP30 ENERGY']}
pv1['SAPT2+(3) TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) ELST ENERGY', 'SAPT2+(3) EXCH ENERGY', 'SAPT2+(3) INDC ENERGY', 'SAPT2+(3) DISP ENERGY']}
pv1['SAPT2+(3)(CCD) ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) ELST ENERGY']}
pv1['SAPT2+(3)(CCD) EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) EXCH ENERGY']}
pv1['SAPT2+(3)(CCD) INDC ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) INDC ENERGY']}
pv1['SAPT2+(3)(CCD) DISP ENERGY'] = {'func': sum, 'args': ['SAPT CCD DISP', 'SAPT DISP30 ENERGY']}
pv1['SAPT2+(3)(CCD) TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(3)(CCD) ELST ENERGY', 'SAPT2+(3)(CCD) EXCH ENERGY', 'SAPT2+(3)(CCD) INDC ENERGY', 'SAPT2+(3)(CCD) DISP ENERGY']}
pv1['SAPT2+(3)DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) ELST ENERGY']}
pv1['SAPT2+(3)DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) EXCH ENERGY']}
pv1['SAPT2+(3)DMP2 INDC ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) INDC ENERGY', 'SAPT MP2(2) ENERGY']}
pv1['SAPT2+(3)DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) DISP ENERGY']}
pv1['SAPT2+(3)DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(3)DMP2 ELST ENERGY', 'SAPT2+(3)DMP2 EXCH ENERGY', 'SAPT2+(3)DMP2 INDC ENERGY', 'SAPT2+(3)DMP2 DISP ENERGY']}
pv1['SAPT2+(3)(CCD)DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) ELST ENERGY']}
pv1['SAPT2+(3)(CCD)DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) EXCH ENERGY']}
pv1['SAPT2+(3)(CCD)DMP2 INDC ENERGY'] = {'func': sum, 'args': ['SAPT2+(3)DMP2 INDC ENERGY']}
pv1['SAPT2+(3)(CCD)DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+(3)(CCD) DISP ENERGY']}
pv1['SAPT2+(3)(CCD)DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(3)(CCD)DMP2 ELST ENERGY', 'SAPT2+(3)(CCD)DMP2 EXCH ENERGY', 'SAPT2+(3)(CCD)DMP2 INDC ENERGY', 'SAPT2+(3)(CCD)DMP2 DISP ENERGY']}
pv1['SAPT2+3 ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY', 'SAPT ELST12,R ENERGY', 'SAPT ELST13,R ENERGY']}
pv1['SAPT2+3 EXCH ENERGY'] = {'func': lambda x: x[1] + x[0] * (x[2] + x[3]), 
                              'args': ['SAPT EXCHSCAL', 'SAPT EXCH10 ENERGY', 'SAPT EXCH11(S^2) ENERGY', 'SAPT EXCH12(S^2) ENERGY']}
pv1['SAPT2+3 INDC ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5] + x[6] + x[0] * x[7],
                              'args': ['SAPT EXCHSCAL', 'SAPT HF(3) ENERGY', 'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY',
                                       'SAPT IND22 ENERGY', 'SAPT EXCH-IND22 ENERGY', 'SAPT IND30,R ENERGY', 'SAPT EXCH-IND30,R ENERGY']}
pv1['SAPT2+3 DISP ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5],
                              'args': ['SAPT EXCHSCAL', 'SAPT MP4 DISP', 'SAPT DISP30 ENERGY', 'SAPT EXCH-DISP30 ENERGY',
                                       'SAPT IND-DISP30 ENERGY', 'SAPT EXCH-IND-DISP30 ENERGY']}
pv1['SAPT2+3 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+3 ELST ENERGY', 'SAPT2+3 EXCH ENERGY', 'SAPT2+3 INDC ENERGY', 'SAPT2+3 DISP ENERGY']}
pv1['SAPT2+3(CCD) ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+3 ELST ENERGY']}
pv1['SAPT2+3(CCD) EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+3 EXCH ENERGY']}
pv1['SAPT2+3(CCD) INDC ENERGY'] = {'func': sum, 'args': ['SAPT2+3 INDC ENERGY']}
pv1['SAPT2+3(CCD) DISP ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5],
                                   'args': ['SAPT EXCHSCAL', 'SAPT CCD DISP', 'SAPT DISP30 ENERGY', 'SAPT EXCH-DISP30 ENERGY',
                                         'SAPT IND-DISP30 ENERGY', 'SAPT EXCH-IND-DISP30 ENERGY']}
pv1['SAPT2+3(CCD) TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+3(CCD) ELST ENERGY', 'SAPT2+3(CCD) EXCH ENERGY', 'SAPT2+3(CCD) INDC ENERGY', 'SAPT2+3(CCD) DISP ENERGY']}
pv1['SAPT2+3DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+3 ELST ENERGY']}
pv1['SAPT2+3DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+3 EXCH ENERGY']}
pv1['SAPT2+3DMP2 INDC ENERGY'] = {'func': sum, 'args': ['SAPT2+3 INDC ENERGY', 'SAPT MP2(3) ENERGY']}
pv1['SAPT2+3DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+3 DISP ENERGY']}
pv1['SAPT2+3DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+3DMP2 ELST ENERGY', 'SAPT2+3DMP2 EXCH ENERGY', 'SAPT2+3DMP2 INDC ENERGY', 'SAPT2+3DMP2 DISP ENERGY']}
pv1['SAPT2+3(CCD)DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+3 ELST ENERGY']}
pv1['SAPT2+3(CCD)DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+3 EXCH ENERGY']}
pv1['SAPT2+3(CCD)DMP2 INDC ENERGY'] = {'func': sum, 'args': ['SAPT2+3DMP2 INDC ENERGY']}
pv1['SAPT2+3(CCD)DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+3(CCD) DISP ENERGY']}
pv1['SAPT2+3(CCD)DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+3(CCD)DMP2 ELST ENERGY', 'SAPT2+3(CCD)DMP2 EXCH ENERGY', 'SAPT2+3(CCD)DMP2 INDC ENERGY', 'SAPT2+3(CCD)DMP2 DISP ENERGY']}
pv1['DFT-SAPT ELST ENERGY'] = {'func': sum, 'args': ['DFT-SAPT ELST10,R ENERGY']}
pv1['DFT-SAPT EXCH ENERGY'] = {'func': sum, 'args': ['DFT-SAPT EXCH10 ENERGY']}
pv1['DFT-SAPT INDC ENERGY'] = {'func': sum, 'args': ['SAPT HF(2) ALPHA=0.0 ENERGY', 
                                                     'DFT-SAPT IND20,R ENERGY', 'DFT-SAPT EXCH-IND20,R ENERGY']}
pv1['DFT-SAPT DISP ENERGY'] = {'func': sum, 'args': ['DFT-SAPT DISP20 ENERGY', 'DFT-SAPT EXCH-DISP20 ENERGY']}
pv1['DFT-SAPT TOTAL ENERGY'] = {'func': sum, 'args': ['DFT-SAPT ELST ENERGY', 'DFT-SAPT EXCH ENERGY', 'DFT-SAPT INDC ENERGY', 'DFT-SAPT DISP ENERGY']}
for pvar, action in pv1.iteritems():
    try:
        if verbose > 0:
            print """building %s %s""" % (pvar, '.' * (50 - len(pvar))),
        df = append_result_of_func_with_funcargs_to_master_DataFrame_atlevel_with_label(
            master=df, label=pvar, atlevel=lvl, func=action['func'], funcargs=action['args'])
        if verbose > 0:
            print """SUCCESS"""
    except KeyError, e:
        if verbose > 0:
            print """FAILED, missing %s""" % (e)

if project == 'merz3':
    qlvl = 'SAPT0'
    qbas = 'jadz'
    print '{}/{} for {} in {}:  RXN, TOT, ELST, EXCH, IND, DISP'.format(qlvl, qbas, len(rxns), dbse)
    for rxn in rxns:
        qtotl = h2kc * df.xs('{} TOTAL ENERGY'.format(qlvl), level='psivar').xs(qbas, level='bstrt').xs(rxn, level='rxn')['CP']['dimer']
        qelst = h2kc * df.xs('{} ELST ENERGY'.format(qlvl), level='psivar').xs(qbas, level='bstrt').xs(rxn, level='rxn')['CP']['dimer']
        qexch = h2kc * df.xs('{} EXCH ENERGY'.format(qlvl), level='psivar').xs(qbas, level='bstrt').xs(rxn, level='rxn')['CP']['dimer']
        qindc = h2kc * df.xs('{} INDC ENERGY'.format(qlvl), level='psivar').xs(qbas, level='bstrt').xs(rxn, level='rxn')['CP']['dimer']
        qdisp = h2kc * df.xs('{} DISP ENERGY'.format(qlvl), level='psivar').xs(qbas, level='bstrt').xs(rxn, level='rxn')['CP']['dimer']
        print """%-23s %8.4f   %8.4f %8.4f %8.4f %8.4f""" % (rxn, qtotl, qelst, qexch, qindc, qdisp)

    qlvl = 'SSAPT0'
    print '{}/{} for {} in {}:  RXN, TOT, ELST, EXCH, IND, DISP'.format(qlvl, qbas, len(rxns), dbse)
    for rxn in rxns:
        qtotl = h2kc * df.xs('{} TOTAL ENERGY'.format(qlvl), level='psivar').xs(qbas, level='bstrt').xs(rxn, level='rxn')['CP']['dimer']
        qelst = h2kc * df.xs('{} ELST ENERGY'.format(qlvl), level='psivar').xs(qbas, level='bstrt').xs(rxn, level='rxn')['CP']['dimer']
        qexch = h2kc * df.xs('{} EXCH ENERGY'.format(qlvl), level='psivar').xs(qbas, level='bstrt').xs(rxn, level='rxn')['CP']['dimer']
        qindc = h2kc * df.xs('{} INDC ENERGY'.format(qlvl), level='psivar').xs(qbas, level='bstrt').xs(rxn, level='rxn')['CP']['dimer']
        qdisp = h2kc * df.xs('{} DISP ENERGY'.format(qlvl), level='psivar').xs(qbas, level='bstrt').xs(rxn, level='rxn')['CP']['dimer']
        print """%-23s %8.4f   %8.4f %8.4f %8.4f %8.4f""" % (rxn, qtotl, qelst, qexch, qindc, qdisp)
    

#for ml in mlist:
#    mdf = df.xs(ml, level='psivar').xs('jadz', level='bstrt') #.xs('BBI-004GLU-063LEU-2', level='rxn')
#    #print ml, '\n'#, categories(mdf, 0)
#    mdf *= h2kc
#    print mdf.head(5)

#mlist = ['SAPT0 TOTAL ENERGY', 'SAPT0 DISP ENERGY', 'SAPT EXCHSCAL', 'SAPT HF(2) ENERGY']
#for ml in mlist:
#    mdf = df.xs(ml, level='psivar').xs('adz', level='bstrt').xs('S22-2', level='rxn')
#    print ml, '\n'#, categories(mdf, 0)
#    #print h2kc * mdf
#    print mdf.head(5)

#print df.xs('B3LYP TOTAL ENERGY', level='psivar').head(10)
#print df.xs('B3LYP-D3 TOTAL ENERGY', level='psivar').head(10)
#print df.xs('B3LYP TOTAL ENERGY', level='psivar').xs('dfhf', level='meta').head(10)
#print df.xs('B3LYP-D3 TOTAL ENERGY', level='psivar').xs('dfhf', level='meta').head(10)

# <<< define extra Method and BasisSet objects >>>

# Note dict key must match object name; i.e., first two strings on line must be the same
#   Note also that methods are all caps, bases all lowercase

#bases['hill1_adtz'] = BasisSet('hill1_adtz', build=[['hillcc_adtz'], ['atz', 'hillcc_adtz']])  # TODO should have None or non-xtpl first element?
#bases['hill2_dtzf12'] = BasisSet('hill2_dtzf12', build=[None, None, ['tzf12', 'hillcc_dtzf12', 'hillt_dtzf12']])
#methods['CCSDTNSAF12'] = Method('CCSDTNSAF12', fullname='CCSD(T)-F12a')

# <<< append to main DataFrame computable basis treatment quantities >>>

lvl = 'bstrt'
pv2 = collections.OrderedDict()
pv2['dtz'] = {'func': xtpl_power, 'args': [3.0, 3, 'tz', 'dz']}
pv2['jadtz'] = {'func': xtpl_power, 'args': [3.0, 3, 'jatz', 'jadz']}
pv2['hadtz'] = {'func': xtpl_power, 'args': [3.0, 3, 'hatz', 'hadz']}
pv2['adtz'] = {'func': xtpl_power, 'args': [3.0, 3, 'atz', 'adz']}
pv2['tqz'] = {'func': xtpl_power, 'args': [3.0, 4, 'qz', 'tz']}
pv2['matqz'] = {'func': xtpl_power, 'args': [3.0, 4, 'maqz', 'matz']}
pv2['jatqz'] = {'func': xtpl_power, 'args': [3.0, 4, 'jaqz', 'jatz']}
pv2['hatqz'] = {'func': xtpl_power, 'args': [3.0, 4, 'haqz', 'hatz']}
pv2['atqz'] = {'func': xtpl_power, 'args': [3.0, 4, 'aqz', 'atz']}
pv2['q5z'] = {'func': xtpl_power, 'args': [3.0, 5, '5z', 'qz']}
pv2['maq5z'] = {'func': xtpl_power, 'args': [3.0, 5, 'ma5z', 'maqz']}
pv2['jaq5z'] = {'func': xtpl_power, 'args': [3.0, 5, 'ja5z', 'jaqz']}
pv2['haq5z'] = {'func': xtpl_power, 'args': [3.0, 5, 'ha5z', 'haqz']}
pv2['aq5z'] = {'func': xtpl_power, 'args': [3.0, 5, 'a5z', 'aqz']}
pv2['56z'] = {'func': xtpl_power, 'args': [3.0, 6, '6z', '5z']}
pv2['ma56z'] = {'func': xtpl_power, 'args': [3.0, 6, 'ma6z', 'ma5z']}
pv2['ja56z'] = {'func': xtpl_power, 'args': [3.0, 6, 'ja6z', 'ja5z']}
pv2['ha56z'] = {'func': xtpl_power, 'args': [3.0, 6, 'ha6z', 'ha5z']}
pv2['a56z'] = {'func': xtpl_power, 'args': [3.0, 6, 'a6z', 'a5z']}
pv2['dtzf12'] = {'func': xtpl_power, 'args': [3.0, 3, 'tzf12', 'dzf12']}
pv2['tqzf12'] = {'func': xtpl_power, 'args': [3.0, 4, 'qzf12', 'tzf12']}
# Hill xtpl for CCSD-F12b from Table X of JCP 131 194105 (2009)
    # TODO should only be applied to CCF12, as per definition (literally, only CCF12B)
pv2['hillcc_adtz'] = {'func': xtpl_power, 'args': [2.483070, 3, 'atz', 'adz']}
pv2['hillcc_atqz'] = {'func': xtpl_power, 'args': [4.255221, 4, 'aqz', 'atz']}
pv2['hillcc_aq5z'] = {'func': xtpl_power, 'args': [4.910269, 5, 'a5z', 'aqz']}
pv2['hillcc_dtzf12'] = {'func': xtpl_power, 'args': [3.144518, 3, 'tzf12', 'dzf12']}
pv2['hillcc_tqzf12'] = {'func': xtpl_power, 'args': [4.595995, 4, 'qzf12', 'tzf12']}
# Hill xtpl for unscaled (T)-F12 from Table XI of JCP 131 194105 (2009)
pv2['hillt_adtz'] = {'func': xtpl_power, 'args': [2.790300, 3, 'atz', 'adz']}
pv2['hillt_dtzf12'] = {'func': xtpl_power, 'args': [2.615472, 3, 'tzf12', 'dzf12']}
for pvar, action in pv2.iteritems():
    try:
        if verbose > 0:
            print """building %s %s""" % (pvar, '.' * (50 - len(pvar))),
        df = append_result_of_func_with_funcargs_to_master_DataFrame_atlevel_with_label(
            master=df, label=pvar, atlevel=lvl, func=action['func'], funcargs=action['args'])
        if verbose > 0:
            print """SUCCESS"""
    except KeyError, e:
        if verbose > 0:
            print """FAILED, missing %s""" % (e)




#     <<< SCS(MI)-MP2 & SCS(MI)-F12 >>>
try:
    df.xs('MP2-F12 CORRELATION ENERGY', level='psivar')
    df.xs('MP2 CORRELATION ENERGY', level='psivar')
    pass  # TODO
except KeyError, e:
    print 'Not handled: SCS(MI)-MP2 OMEGA', e
else:
    mi_os = {'atz': 0.17,  'adtz': 0.29,  'aqz': 0.31,  'atqz': 0.40,
            'hatz': 0.17, 'hadtz': 0.29, 'haqz': 0.31, 'hatqz': 0.40,
            'jatz': 0.17, 'jadtz': 0.29, 'jaqz': 0.31, 'jatqz': 0.40,
            'matz': 0.17,                'maqz': 0.31, 'matqz': 0.40,
                                         'aaqz': 0.31,
              'tz': 0.17,   'dtz': 0.29,   'qz': 0.31,   'tqz': 0.40,
           'tzf12': 0.17}
    mi_ss = {'atz': 1.75,  'adtz': 1.46,  'aqz': 1.46,  'atqz': 1.29,
            'hatz': 1.75, 'hadtz': 1.46, 'haqz': 1.46, 'hatqz': 1.29,
            'jatz': 1.75, 'jadtz': 1.46, 'jaqz': 1.46, 'jatqz': 1.29,
            'matz': 1.75,                'maqz': 1.46, 'matqz': 1.29,
                                         'aaqz': 1.46,
              'tz': 1.75,   'dtz': 1.46,   'qz': 1.46,   'tqz': 1.29,
           'tzf12': 1.75}
    mi_os = collections.defaultdict(lambda: np.nan, mi_os)
    mi_ss = collections.defaultdict(lambda: np.nan, mi_ss)
    merge = {}
    merge_os = {}
    merge_ss = {}
    for bstrt in categories(df.xs('MP2 CORRELATION ENERGY', level='psivar'), 0):
        temp_os = df.xs('MP2 CORRELATION ENERGY', level='psivar').xs(bstrt, level='bstrt')
        temp_ss = df.xs('MP2 CORRELATION ENERGY', level='psivar').xs(bstrt, level='bstrt')
        temp_os[:][:] = mi_os[bstrt]
        temp_ss[:][:] = mi_ss[bstrt]
        merge_os[bstrt] = temp_os
        merge_ss[bstrt] = temp_ss
    merge['SCS(MI)-MP2 SCS-OS'] = pd.concat(merge_os)
    merge['SCS(MI)-MP2 SCS-SS'] = pd.concat(merge_ss)
    merge_os = {}
    merge_ss = {}
    for bstrt in categories(df.xs('MP2-F12 CORRELATION ENERGY', level='psivar'), 0):
        temp_os = df.xs('MP2-F12 CORRELATION ENERGY', level='psivar').xs(bstrt, level='bstrt')
        temp_ss = df.xs('MP2-F12 CORRELATION ENERGY', level='psivar').xs(bstrt, level='bstrt')
        temp_os[:][:] = mi_os[bstrt]
        temp_ss[:][:] = mi_ss[bstrt]
        merge_os[bstrt] = temp_os
        merge_ss[bstrt] = temp_ss
    merge['SCS(MI)-MP2-F12 SCS-OS'] = pd.concat(merge_os)
    merge['SCS(MI)-MP2-F12 SCS-SS'] = pd.concat(merge_ss)
    df_mi = pd.concat(merge)
    df_mi = df_mi.reorder_levels([1, 0, 2, 3])
    df_mi.index.names = ['bstrt', 'psivar', 'meta', 'rxn']
    df_mi = df_mi.dropna(how='all')
    if len(df_mi) > 0:
        df = df.append(df_mi, verify_integrity=True)

lvl = 'psivar'
pv3 = collections.OrderedDict()
pv3['SCS(MI)-MP2 CORRELATION ENERGY'] = {'func': spin_component_scaling, 'args': ['SCS(MI)-MP2 SCS-OS', 'SCS(MI)-MP2 SCS-SS', 'MP2 CORRELATION ENERGY', 'MP2 SAME-SPIN CORRELATION ENERGY']}
pv3['SCS(MI)-MP2 TOTAL ENERGY'] = {'func': sum, 'args': ['HF TOTAL ENERGY', 'SCS(MI)-MP2 CORRELATION ENERGY']}
pv3['SCS(MI)-MP2-F12 CORRELATION ENERGY'] = {'func': spin_component_scaling, 'args': ['SCS(MI)-MP2-F12 SCS-OS', 'SCS(MI)-MP2-F12 SCS-SS', 'MP2-F12 CORRELATION ENERGY', 'MP2-F12 SAME-SPIN CORRELATION ENERGY']}
pv3['SCS(MI)-MP2-F12 TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'SCS(MI)-MP2-F12 CORRELATION ENERGY']}
for pvar, action in pv3.iteritems():
    try:
        if verbose > 0:
            print """building %s %s""" % (pvar, '.' * (50 - len(pvar))),
        df = append_result_of_func_with_funcargs_to_master_DataFrame_atlevel_with_label(
            master=df, label=pvar, atlevel=lvl, func=action['func'], funcargs=action['args'])
        if verbose > 0:
            print """SUCCESS"""
    except KeyError, e:
        if verbose > 0:
            print """FAILED, missing %s""" % (e)


# <<< define simple functions codifying cbs() piecing >>>

# TODO move these to be functions in Method or BasisSet

# min and max pieces for which basis treatment intent is clearly understood
# with the understanding that Method should be defining all those pieces that have uniform basis
#           desi max_mtd    min_bas    max_bas
#HF/adz     1    1          1           inf (1)
#MP2/adz    1    2          1           inf (1)
#CC/adz     1    3          1           inf (1)
#HF/adtz    1    1          1           2
#MP2/adtz   2    2          1           2
#CC/adtz    2    3          1           2
#CC/atqzadz 3    3          2           3
#   stages = min(max_mtd, max_bas)


def compute_max_bas(bas):
    return len(bases[bas].build)


def generic_bas(Nstage, bas):
    return bases[bas].build[Nstage - 1]  # 1-indexed stage counting to 0-indexed array storage


def compute_max_mtd(stub):
    if stub in ['HFCABS', 'HF', 'SCF']:
        return 1
    elif 'MP2' in stub and stub not in ['MP25', 'MP2C', 'MP2CF12']:
        return 2
    else:
        return 3


def generic_mtd(Nstage, stub):
    if Nstage == 1:
        return ['%s TOTAL ENERGY' % (stub)]
    elif Nstage == 2:
        return ['%s TOTAL ENERGY' % ('HF-CABS' if 'F12' in stub else 'SCF'), 
                '%s CORRELATION ENERGY' % (stub)]
    elif Nstage == 3:
        #return ['%s TOTAL ENERGY' % ('HF-CABS' if 'F12' in stub else 'SCF'), '%s CORRELATION ENERGY' % (stub), '%s CC CORRECTION ENERGY' % (stub)]
        return ['%s TOTAL ENERGY' % ('HF-CABS' if 'F12' in stub else 'SCF'), 
                '%s CORRELATION ENERGY' % ('MP2-F12' if 'F12' in stub else 'MP2'),
                '%s CC CORRECTION ENERGY' % (stub)]


def build_from_lists(mtdlist, baslist, optlist=None):
    if optlist is None:
        optlist = [''] * len(mtdlist)
    return ie(sum([df.loc[bas].loc[pcs].loc[opt] for pcs, bas, opt in zip(mtdlist, baslist, optlist)]))


def build(method, option, cpmode, basis):
    #print method, option, cpmode, basis
    Nstage = min(compute_max_mtd(method), compute_max_bas(basis))
    baslist = generic_bas(Nstage, basis)
    mtdlist = generic_mtd(Nstage, methods[method].fullname.upper())
    optlist = ['-'.join([opt for opt in option.split('-') if (opt == '' or pcs in optclue2psivar[opt])]) for pcs in mtdlist]
    # not sure about how this will hold up to multistage
    #   purpose is to prevent duplicate IE with irrelevant option labels surviving
    for alw in option.split('-'):
        found = False
        for piece in optlist:
            if alw in piece:
                found = True
        if not found:
            return 'Overly'
    func = {'CP': ie2, 'unCP': ie_uncp2, 'ave': ie_ave, 'default': default}[cpmode]
    if baslist is None:
        raise KeyError  # TODO a more specific message that mtd/bas don't mix wouldn't hurt
        #print '\n <<<', methods[method].fullname, '/', bases[basis].fullname, '>>>'
        #print 'stages:', 'M:', compute_max_mtd(method), 'B:', compute_max_bas(basis), 'U:', Nstage
    if method in [] and basis in []:
    #if method in ['MP2'] and basis in ['aqz']:
    #if method in ['MP2C', 'MP2CF12'] and basis in ['atqzadz', 'adz']:
    #if method in ['B3LYPD3'] and basis in ['def2qzvp']:
        print '\n', method, option, cpmode, basis
        print 'pcss:', mtdlist
        print 'bass:', baslist
        print 'opts:', optlist
        for pcs, bas, opt in zip(mtdlist, baslist, optlist):
            print df.loc[bas].loc[pcs].loc[opt].loc['ACONF-15'] #'NBC1-BzBz_S-5.0'] #'S22-2'] #'A24-1'] #'BBI-150LYS-158LEU-2'] #'S22-2']
    return func(sum([df.loc[bas].loc[pcs].loc[opt] for pcs, bas, opt in zip(mtdlist, baslist, optlist)]))

# <<< assemble all model chemistries into columns of new DataFrame >>>

rxns = ['%s-%s' % (dbse, rxn) for rxn in dbobj.dbdict[dbobj.dbse].hrxn.keys()]
mine = pd.DataFrame({}, index=rxns)
mine.index.names = ['rxn']

if project == 'dft':
    mtds = ['B3LYP', 'B3LYPD2', 'B3LYPD3', 'B2PLYP', 'B2PLYPD2', 'B2PLYPD3',
            'B970', 'BP86', 'B97', 'WB97XD', 'M052X', 'M062X', 'PBE', 'PBE0', 'XYG3',
            'B970D2', 'BP86D2', 'BP86D3', 'M052XD3', 'M062XD3', 'PBED2', 'PBED3',
            'PBE0D2', 'PBE0D3', 'B97D2', 'B97D3', 'B3LYPXDM']
    bass = ['adz', 'atz', 'dadz', 'datz', 'dz', 'tz', 'hadz', 'hatz',
            '6311pg_3df_2p_', '6311ppg_3df_2p_']
    opts = ['']
    cpmd = ['unCP', 'CP']

elif project == 'parenq':
    mtds = ['HF', 'CCSDT', 'CCSDTQ', 'CCSD', 'CCSDFULLT', 'MP2']
    bass = ['adz', 'atz', 'hadz', 'jadz', 'atz', 'adtz', 'aq5z']
    opts = ['', 'full', 'fno1e5', 'fno1e4', 'fno5e5', 'fno1e6', 'fno1e4-mrcc', 'fno1e3']
    cpmd = ['CP']

elif project == 'f12dilabio':
    mtds = [mtd for mtd in methods if ((mtd == 'HFCABS' or 'F12' in mtd) and ('SC' not in mtd and 'MP2C' not in mtd and 'DWMP2' not in mtd))]
    bass = ['adz', 'atz', 'aqz', 'a5z', 'dzf12', 'tzf12', 'qzf12',
            'adtz', 'atqz', 'aq5z', 'dtzf12', 'tqzf12',
            'hill1_adtz', 'hill1_atqz', 'hill1_aq5z', 'hill1_dtzf12', 'hill1_tqzf12']
    opts = ['']
    cpmd = ['CP']

elif project == 'dilabio':
    mtds = ['HF', 'MP2', 'CCSD', 'CCSDT']
    bass = ['adz', 'atz', 'aqz', 'a5z', 'a6z', 'adtz', 'atqz', 'aq5z', 'a56z',
            'atqzadz', 'atqzatz', 'aq5zadz', 'aq5zatz']
    opts = ['']
    cpmd = ['CP', 'unCP', 'ave']

elif project == 'dhdft':
    mtds = ['B3LYP', 'B3LYPD3', 'B2PLYP', 'B2PLYPD3', 'B97D3', 'M052X', 'M062X',
            'DSDPBEP86', 'VV10', 'LCVV10', 'DSDPBEP86', 'PBE', 'PBED3', 'PBE0', 'PBE0D3',
            'M08HX', 'M08SO', 'M11', 'M11L', 'PBE02', 'DLDFD', 'WB97X2', 'WB97XD',
            'M062XD3']
    bass = ['adz', 'atz']
    opts = ['', 'nfc']
    cpmd = ['CP', 'unCP']
    
elif project == 'dhdft2':
    mtds = ['DSDPBEP86', 'DSDPBEP86D2', 'DSDPBEP86D3BJ']
    bass = ['adz', 'atz']
    opts = ['', 'nfc']
    cpmd = ['CP', 'unCP']
    
elif project == 'pt2':
    mtds = ['HF', 'MP2', 'SCSMP2', 'SCSNMP2', 'SCSMIMP2', 'DWMP2', 'MP2C',
            'MP25', 'MP3',
            'HFCABS', 'MP2F12', 'SCSMP2F12', 'SCSNMP2F12', 'SCSMIMP2F12', 'DWMP2F12', 'MP2CF12',
            'CCSDAF12', 'SCSCCSDAF12', 'SCMICCSDAF12', 'CCSDTAF12', 'CCSDTSAF12', 'CCSDTNSAF12',
            'CCSDBF12', 'SCSCCSDBF12', 'SCMICCSDBF12', 'CCSDTBF12', 'CCSDTSBF12', 'CCSDTNSBF12', 
            'DWCCSDTF12', 'DWCCSDTSF12', 'DWCCSDTNSF12',        
            'SAPT0', 'SAPT0S', 'SAPTSCS', 'SAPTDFT', 'SAPT2',
            'SAPT2P', 'SAPT2PC', 'SAPT2PM', 'SAPT2PCM',
            'SAPT3', 'SAPT3C', 'SAPT3M', 'SAPT3CM',
            'SAPT3F', 'SAPT3FC', 'SAPT3FM', 'SAPT3FCM',
            ]
    bass = ['adz', 'atz', 'aqz', 'a5z',
            'adtz', 'atqz', 'aq5z',
            'atzadz', 'adtzadz', 'atqzadz', 'atqzatz', 'atzadtz', 'atqzadtz',
            'dz', 'jadz', 'hadz',
            'tz', 'matz', 'jatz', 'hatz',
            'qz', 'aaqz', 'maqz', 'jaqz', 'haqz',
            'dtz', 'jadtz', 'hadtz',
            'tqz', 'matqz', 'jatqz', 'hatqz',
            'atzdz', 'adtzdz', 'atqzdz', 'atqztz', 'atzdtz', 'atqzdtz',
            'atzjadz', 'adtzjadz', 'atqzjadz', 'atqzjatz', 'atzjadtz', 'atqzjadtz',
            'atzhadz', 'adtzhadz', 'atqzhadz', 'atqzhatz', 'atzhadtz', 'atqzhadtz',
            'dzf12', 'tzf12']
    opts = ['', 'dfhf', 'dfmp', 'dfhf-dfmp']
    cpmd = ['CP']
    
elif project == 'saptone':
    mtds = ['SAPT0', 'SAPT0S', 'SAPTSCS', 'SAPTDFT', 'SAPT2',
            'SAPT2P', 'SAPT2PC', 'SAPT2PM', 'SAPT2PCM',
            'SAPT3', 'SAPT3C', 'SAPT3M', 'SAPT3CM',
            'SAPT3F', 'SAPT3FC', 'SAPT3FM', 'SAPT3FCM']
    bass = ['adz', 'atz', 'aqz', 'a5z',
            'dz', 'jadz', 'hadz',
            'tz', 'matz', 'jatz', 'hatz',
            'qz', 'aaqz', 'maqz', 'jaqz', 'haqz']
    opts = ['', 'dfmp']
    cpmd = ['CP']

elif project == 'merz3':
    #mtds = ['MP2', 'SCSMP2', 'SCSNMP2', 'SCSMIMP2', 'DWMP2', 'MP2C',
    #        'MP2CF12', 'DWCCSDTF12']
    #bass = ['adz', 'atz', 'aqz', 'addz', 'adtz', 'atqz', 'qz', 'atqzadz']
    #opts = ['', 'dfhf-dfmp']
    mtds = ['SAPT0', 'SAPT0S']
    bass = ['jadz']
    opts = ['']  # should SAPT0 be coming through with a dfmp label?
    cpmd = ['CP']
    
elif project == 'aep':
    mtds = ['SCSMICCSD']
    bass = ['adz', 'atz', 'adtz', 'atqzadz', 'atqzatz']
    opts = ['', 'dfhf-dfmp']
    cpmd = ['CP']

elif project == 'sflow':
    mtds = ['MP2']
    bass = ['aqz']
    opts = ['dfhf-dfmp-dsrgs0p1', 'dfhf-dfmp-dsrgs0p5', 'dfhf-dfmp-dsrgs1p0']
    cpmd = ['CP']

elif project == 'dfit':
    mtds = ['B3LYP', 'B97', 'BLYP', 'BP86', 'PBE', 'PBE0', 'WPBE', 'B2PLYP',
            'B3LYPD3', 'B97D3', 'BLYPD3', 'BP86D3', 'PBED3', 'PBE0D3', 'B2PLYPD3',]
    bass = ['def2qzvp']
    opts = ['', 'dfhf', 'dfhf-dfmp']  # why B2PLYP with no or dfhf label; why B3LYP w dfhf-dfmp label?
    cpmd = ['CP', 'unCP']

elif project == 'nbcref':
    mtds = ['CCSDT']
    bass = ['atqzhatz', 'atqzatz', 'atqzadz', 'atqz']
    opts = ['', 'dfhf-dfmp']
    cpmd = ['CP']

elif project == 'curveref':
    mtds = ['CCSDT', 'DWCCSDTF12',
            'SCMICCSDAF12', 'CCSDTAF12', 'SCMICCSDBF12', 'CCSDTBF12', 'SCMICCSDCF12', 'CCSDTCF12']
    bass = ['atz', 'aqz', 'a5z', 'atqz', 'aq5z']
    opts = ['']
    cpmd = ['CP']

for cpm in cpmd:
    for mtd in mtds:
        for bas in bass:
            for opt in opts:
                mc = '-'.join([mtd, cpm, bas]) if opt == '' else '-'.join([mtd, opt, cpm, bas])
                try:
                    tmp = build(mtd, opt, cpm, bas)
                    if verbose > 1:
                        print """Built model chemistry %s""" % (mc)
                except KeyError, e:
                    if verbose > 1:
                        print """Error building rxn mc: empty '%s' b/c no %s""" % (mc, e)
                else:
                    if isinstance(tmp, basestring) and tmp == 'Overly':
                        if verbose > 1:
                            print """Overly optioned rxn mc: {} for {}""".format(opt, mc)
                    elif tmp.empty or tmp.dropna(how='all').empty:  # invalid for canopy python
                        if verbose > 1:
                            print """Empty rxn mc: empty '%s'""" % (mc)
                    else:
                        mine[mc] = tmp


def threadtheframe(modelchem, xlimit=4.0):
    dbdat = []
    for rxn in dbobj.hrxn.keys():
        data = []
        for mc in modelchem:
            value = h2kc * mine[mc]['%s-%s' % (dbse, rxn)]
            data.append(None if pd.isnull(value) else value)
        dbdat.append({'sys': str(rxn), 'color': dbobj.dbdict[dbobj.dbse].hrxn[rxn].color, 'data': data})
    qcdb.mpl.thread(dbdat, modelchem, color='sapt', xlimit=xlimit)

if project == 'f12dilabio':
    mine['CCSDTNSBF12-CP-hill2_adtz'] = build_from_lists(['HF-CABS TOTAL ENERGY', 'CCSD-F12B CORRELATION ENERGY', '(T)-F12AB CORRECTION ENERGY'], ['atz', 'hillcc_adtz', 'hillt_adtz'])

if project == 'parenq':
    mine['DELTQ-full-CP-hadz'] = mine['CCSDTQ-full-CP-hadz'] - mine['CCSDT-full-CP-hadz']
    mine['DELTQ-full-CP-jadz'] = mine['CCSDTQ-full-CP-jadz'] - mine['CCSDT-full-CP-jadz']

    mine['DELTQ-full-CP-adz'] = mine['CCSDTQ-full-CP-adz'] - mine['CCSDT-full-CP-adz']
    mine['DELTQ-full-CP-atz'] = mine['CCSDTQ-full-CP-atz'] - mine['CCSDT-full-CP-atz']
    mine['DELTQ-full-CP-adtz'] = mine['CCSDTQ-full-CP-adtz'] - mine['CCSDT-full-CP-adtz']
    
    mine['DELTQ-fno1e4-CP-adz'] = mine['CCSDTQ-fno1e4-CP-adz'] - mine['CCSDT-fno1e4-CP-adz']
    mine['DELTQ-fno1e4-CP-atz'] = mine['CCSDTQ-fno1e4-CP-atz'] - mine['CCSDT-fno1e4-CP-atz']
    mine['DELTQ-fno1e4-CP-adtz'] = mine['CCSDTQ-fno1e4-CP-adtz'] - mine['CCSDT-fno1e4-CP-adtz']

    mine['DELTQ-fno1e5-CP-adz'] = mine['CCSDTQ-fno1e5-CP-adz'] - mine['CCSDT-fno1e5-CP-adz']
    mine['DELTQ-fno1e5-CP-atz'] = mine['CCSDTQ-fno1e5-CP-atz'] - mine['CCSDT-fno1e5-CP-atz']
    mine['DELTQ-fno1e5-CP-adtz'] = mine['CCSDTQ-fno1e5-CP-adtz'] - mine['CCSDT-fno1e5-CP-adtz']

    #mine['DELTQ-fno1e3-CP-adz'] = mine['CCSDTQ-fno1e3-CP-adz'] - mine['CCSDT-fno1e3-CP-adz']
    mine['DELTQ-fno1e6-CP-adz'] = mine['CCSDTQ-fno1e6-CP-adz'] - mine['CCSDT-fno1e6-CP-adz']
    mine['DELTQ-fno5e5-CP-adz'] = mine['CCSDTQ-fno5e5-CP-adz'] - mine['CCSDT-fno5e5-CP-adz']


    mine['DEL2T-full-CP-hadz'] = mine['CCSDT-full-CP-hadz'] - mine['MP2-full-CP-hadz']
    mine['DEL2T-full-CP-jadz'] = mine['CCSDT-full-CP-jadz'] - mine['MP2-full-CP-jadz']

    mine['DEL2T-full-CP-adz'] = mine['CCSDT-full-CP-adz'] - mine['MP2-full-CP-adz']
    mine['DEL2T-full-CP-atz'] = mine['CCSDT-full-CP-atz'] - mine['MP2-full-CP-atz']
    mine['DEL2T-full-CP-adtz'] = mine['CCSDT-full-CP-adtz'] - mine['MP2-full-CP-adtz']
    
    mine['DEL2T-fno1e4-CP-adz'] = mine['CCSDT-fno1e4-CP-adz'] - mine['MP2-fno1e4-CP-adz']
    mine['DEL2T-fno1e4-CP-atz'] = mine['CCSDT-fno1e4-CP-atz'] - mine['MP2-fno1e4-CP-atz']
    mine['DEL2T-fno1e4-CP-adtz'] = mine['CCSDT-fno1e4-CP-adtz'] - mine['MP2-fno1e4-CP-adtz']

    mine['DEL2T-fno1e5-CP-adz'] = mine['CCSDT-fno1e5-CP-adz'] - mine['MP2-fno1e5-CP-adz']
    mine['DEL2T-fno1e5-CP-atz'] = mine['CCSDT-fno1e5-CP-atz'] - mine['MP2-fno1e5-CP-atz']
    mine['DEL2T-fno1e5-CP-adtz'] = mine['CCSDT-fno1e5-CP-adtz'] - mine['MP2-fno1e5-CP-adtz']

    #mine['DEL2T-fno1e3-CP-adz'] = mine['CCSDT-fno1e3-CP-adz'] - mine['MP2-fno1e3-CP-adz']
    mine['DEL2T-fno1e6-CP-adz'] = mine['CCSDT-fno1e6-CP-adz'] - mine['MP2-fno1e6-CP-adz']
    mine['DEL2T-fno5e5-CP-adz'] = mine['CCSDT-fno5e5-CP-adz'] - mine['MP2-fno5e5-CP-adz']
   
if False:
#if project == 'parenq':
    # TODO these aren't actually interaction energies!
    mine['CCSDTQ-corl-CP-adz'] = build_from_lists(['CCSDT(Q) CORRELATION ENERGY'], ['adz'], ['fno1e4'])
    mine['CCSDT-corl-CP-adz'] = build_from_lists(['CCSD(T) CORRELATION ENERGY'], ['adz'], ['fno1e4'])
    mine['CCSD-corl-CP-adz'] = build_from_lists(['CCSD CORRELATION ENERGY'], ['adz'], ['fno1e4-mrcc'])
    mine['MP2-corl-CP-adz'] = build_from_lists(['MP2 CORRELATION ENERGY'], ['adz'])
    mine['CCSDTQ-corl-CP-atz'] = build_from_lists(['CCSDT(Q) CORRELATION ENERGY'], ['atz'], ['fno1e4'])
    mine['CCSDT-corl-CP-atz'] = build_from_lists(['CCSD(T) CORRELATION ENERGY'], ['atz'], ['fno1e4'])
    mine['CCSD-corl-CP-atz'] = build_from_lists(['CCSD CORRELATION ENERGY'], ['atz'], ['fno1e4'])
    mine['MP2-corl-CP-atz'] = build_from_lists(['MP2 CORRELATION ENERGY'], ['atz'])

    mine['convdel-adz'] = mine['CCSDT-corl-CP-adz'] - mine['MP2-corl-CP-adz']
    mine['altdel-adz'] = mine['CCSDTQ-corl-CP-adz'] - mine['CCSD-corl-CP-adz']
    minelist = ['convdel-adz', 'altdel-adz', 'CCSDTQ-corl-CP-adz', 'CCSDT-corl-CP-adz', 'CCSD-corl-CP-adz', 'MP2-corl-CP-adz']
    mine['convdel-atz'] = mine['CCSDT-corl-CP-atz'] - mine['MP2-corl-CP-atz']
    mine['altdel-atz'] = mine['CCSDTQ-corl-CP-atz'] - mine['CCSD-corl-CP-atz']
    minelistatz = ['convdel-atz', 'altdel-atz', 'CCSDTQ-corl-CP-atz', 'CCSDT-corl-CP-atz', 'CCSD-corl-CP-atz', 'MP2-corl-CP-atz']

    import qcdb.mpl
    import matplotlib.pyplot as plt
    #threadtheframe(minelist)
    #threadtheframe(minelistatz)
    #threadtheframe(['CCSDTQ-fno1e4-CP-adz', 'CCSDTQ-fno1e4-CP-atz', 'CCSDTQ-fno1e4-CP-adtz', 'CCSDTQ-fno1e5-CP-adtz', 'CCSDTQ-fno1e5-CP-atz', 'CCSDTQ-fno1e5-CP-adz', 'CCSDT-CP-aq5z'], xlimit=7.0)

if verbose > 1:
    print 'Reaction Energies', mine.columns

# <<< test cases >>>

from qcdb.pandas_test import test_pandas
test_pandas(h2kc, project, mine)

# <<< write qcdb data loader >>>

f1 = open('%s/%s_%s.py' % (homewrite, dbse, project), 'w')
print 'Writing to %s/%s_%s.py ...' % (homewrite, dbse, project)
f1.write('\ndef load_%s(dbinstance):\n\n' % (project))

for mc in mine.columns:
    lmc = mc.split('-')  # TODO could be done better
    method = lmc[0]
    bsse = '_'.join(lmc[1:-1])
    basis = lmc[-1]
    for rxn in dbobj.dbdict[dbobj.dbse].hrxn.keys():
        value = h2kc * mine[mc]['%s-%s' % (dbse, rxn)]
        if pd.isnull(value):
            if rxn == 2:
                print 'dropping', mc
            pass
        else:
            f1.write("""    dbinstance.add_ReactionDatum(dbse='%s', rxn=%s, method='%s', mode='%s', basis='%s', value=%.4f)\n""" %
                (dbse, repr(rxn), method, bsse, basis, value))

f1.close()

# <<< write hdf5 >>>

with pd.get_store('%s/%s_%s.h5' % (homewrite, dbse, project)) as handle:
    handle['pdie'] = mine * h2kc
    print 'hdf5', handle



# <<< collecting section >>>

#print 'indexval before', sorted(set([tup[1] for tup in df.index.values]))

# TODO write test case missing checker into the build loop

# if basis treatment isn't designed for SCF, i.e. adtz where Helgaker xtpl is meant for corl, should first item of build be None?

# form indices for mine from master df
#rxns = []
#[rxns.append(tup[2]) for tup in df.index.values if not rxns.count(tup[2])]

#print 'scf', df.loc['tzf12'].loc['HF-CABS TOTAL ENERGY']
#print 'cc', df.loc['HillCC_dtzf12'].loc['CCSD-F12B CORRELATION ENERGY']
#print '(t)', df.loc['HillT_dtzf12'].loc['(T)-F12AB CORRECTION ENERGY']
#print 'search', df.xs('(T)-F12AB CORRECTION ENERGY', level='psivar')
#print mine['CCSDTNSBF12-CP-Hill2_dtzf12']
#print mine.columns

#print 'TTTTT'
#subject = df.xs('DW-MP2-F12 CORRELATION ENERGY', level='psivar')
##subject.xs('adz', level='bstrt').xs('', level='meta')['CP'][:] = np.nan
#print subject
#print subject.head(5)
#print subject.tail(5)
#print 'len', len(subject), categories(subject, 0), categories(subject, 1)
#subject = subject.dropna(axis=0, how='all')
#print 'TTfrop'
#print subject
#print subject.head(5)
#print subject.tail(5)
#print 'len', len(subject), categories(subject, 0), categories(subject, 1)

#mlist = ['HF-CABS TOTAL ENERGY', #'SCS(MI)-MP2-F12 SCS-OS', 'SCS(MI)-MP2-F12 SCS-SS', 'MP2-F12 CORRELATION ENERGY', 'MP2-F12 SAME-SPIN CORRELATION ENERGY',
#    'SCS(MI)-MP2-F12 CORRELATION ENERGY', 'SCS(MI)-MP2-F12 TOTAL ENERGY']
#for ml in mlist:
#    mdf = df.xs(ml, level='psivar').xs('atz', level='bstrt')
#    print ml, '\n', categories(mdf, 0)
#    print mdf
#    print mdf.head(50)
