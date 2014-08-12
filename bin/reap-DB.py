import os
import sys
import glob
import math
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


# <<< read usemefiles and convert to giant DataFrame >>>

#path = r"""C:\Users\Owner\Documents\f12dilabiousemefiles\usemefiles"""

dbse = 'A24'
project = 'parenq'
path = r"""/Users/loriab/linux/qcdb/data/parenqusemefiles/"""

#dbse = 'S22'
#project = 'dft'
#path = r"""/Users/loriab/linux/qcdb/data/dftusemefiles/"""

#dbse = 'A24'
#project = 'f12dilabio'
#path = r"""/Users/loriab/linux/qcdb/data/f12dilabiousemefiles/"""

#dbse = 'S22'
#project = 'dhdft'
#path = r"""/Users/loriab/linux/qcdb/data/dhdftusemefiles/"""

#dbse = 'A24'
#project = 'dilabio'
#path = r"""/Users/loriab/linux/qcdb/data/dilabiousemefiles/"""

dbobj = qcdb.Database(dbse)
dbse = dbobj.dbse
rxns = ['%s-%s' % (dbse, rxn) for rxn in dbobj.hrxn.keys()]
names = ['rxn', 'dimer', 'monoA', 'monoB']
h2kc = qcdb.psi_hartree2kcalmol
    
rawdata = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(dict)))
usemeglob = glob.glob('%s/%s*useme*' % (path, dbse))
if len(usemeglob) == 0:
    raise ValidationError("""No %s usemefiles in %s.""" % (dbse, path))
for useme in usemeglob:
    spl = os.path.basename(useme).split('.')
    #dbse = spl[0].split('-')[0]
    #ocalc = spl[0].split('-')[1]
    basis = spl[0].split('-')[-1]
    piece = '.'.join(spl[1:])
    optns = spl[0].split('-')[2:-1]
    cpmode = 'unCP'
    optns2 = []
    for opt in sorted(optns):
        if opt in ['CP', 'unCP']:
            cpmode = opt
        else:
            try:
                if useme2psivar[piece] in optclue2psivar[opt]:
                    optns2.append(opt)
            except KeyError as e:
                print 'Error: option %s needs adding to to optclue2psivar is psivarrosetta.py: %s' % (opt, e)
                sys.exit(1)
    optns = '-'.join(optns2)
    #print useme, basis, piece, optns, useme2psivar[piece]
    try:
        print 'Reading useme: %6s %50s %15s %5s' % (basis, useme2psivar[piece], optns, cpmode)
    except KeyError:
        print 'Error: useme %s needs adding to useme2psivar in psivarrosetta.py' % (piece)

    tmp = pd.read_csv('%s' % (useme), index_col=0, sep='\s+', comment='#', na_values='None', names=names)
    if piece.endswith('usemedash'):
        rawdata[basis][useme2psivar[piece]][optns]['unCP'] = tmp.dropna(how='all')
        rawdata[basis][useme2psivar[piece]][optns]['CP'] = tmp.dropna(how='all')
    else:
        rawdata[basis][useme2psivar[piece]][optns][cpmode] = tmp.dropna(how='all')
    #print tmp.head(4)

baszip = {}
for baskey, basval in rawdata.iteritems():
    pvzip = {}
    for pvkey, pvval in basval.iteritems():
        metazip = {}
        for metakey, metaval, in pvval.iteritems():
            print baskey, pvkey, metakey
            metazip[metakey] = pd.concat(metaval, axis=1)  # merge CP/unCP
        pvzip[pvkey] = pd.concat(metazip)
    baszip[baskey] = pd.concat(pvzip)
df = pd.concat(baszip)
df.index.names = ['bstrt', 'psivar', 'meta', 'rxn']


# <<< define utility functions >>>

def ie(rgts):
    cpmode = 'CP'
    return rgts[cpmode]['dimer'] - rgts[cpmode]['monoA'] - rgts[cpmode]['monoB']

def ie_uncp(rgts):
    cpmode = 'unCP'
    return rgts[cpmode]['dimer'] - rgts[cpmode]['monoA'] - rgts[cpmode]['monoB']

def ie_ave(rgts):
    return 0.5 * (rgts['CP']['dimer'] - rgts['CP']['monoA'] - rgts['CP']['monoB'] +
                  rgts['unCP']['dimer'] - rgts['unCP']['monoA'] - rgts['unCP']['monoB'])

def categories(df, lvl):
    return sorted(set([tup[lvl] for tup in df.index.values]))

def append_result_of_func_with_funcargs_to_master_DataFrame_atlevel_with_label(master, label, atlevel, func, funcargs):
    odr = {'psivar': [1, 3, 0, 2],
           'bstrt': [3, 1, 0, 2]}        
    data_rich_args = [master.xs(pv, level=atlevel) if isinstance(pv, basestring) else pv for pv in funcargs]
    multiopt = []
    for item in data_rich_args:
        try:
            multiopt.append(set([tup[1] for tup in item.index.values]))
        except AttributeError:
            pass
    #print '    ... options combinations:', list(itertools.product(*multiopt))
    optzip = {}
    for multiopt in itertools.product(*multiopt):
        imultiopt = iter(multiopt)
        optlabel = '-'.join(set(sorted([opt for opt in multiopt if opt != ''])))
        data_rich_args_2 = [arg.xs(imultiopt.next(), level='meta') if isinstance(arg, pd.DataFrame) else arg for arg in data_rich_args]
        optzip[optlabel] = func(data_rich_args_2)
    temp = pd.concat(optzip)
    temp[atlevel] = label
    temp.set_index(atlevel, append=True, inplace=True)
    #print atlevel, 'QQpre', temp.index.values[0]
    temp = temp.reorder_levels(odr[atlevel])
    #print atlevel, 'WWpst', temp.index.values[0]
    return master.append(temp, verify_integrity=True)


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
        
# <<< append to main DataFrame basic psivar equalities not explicit to useme structure >>>

lvl = 'psivar'
pv0 = collections.OrderedDict()
pv0['MP2-F12 TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv0['CCSD-F12A TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD-F12A CORRELATION ENERGY']}
pv0['CCSD-F12B TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD-F12B CORRELATION ENERGY']}
pv0['CCSD-F12C TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD-F12C CORRELATION ENERGY']}
pv0['CCSD(T**)-F12A TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T**)-F12A CORRELATION ENERGY']}
pv0['CCSD(T**)-F12B TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T**)-F12B CORRELATION ENERGY']}
pv0['CCSD(T**)-F12C TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T**)-F12C CORRELATION ENERGY']}
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
pv0['M05-2X TOTAL ENERGY'] = {'func': sum, 'args': ['M05-2X FUNCTIONAL TOTAL ENERGY']}
pv0['M06-2X TOTAL ENERGY'] = {'func': sum, 'args': ['M06-2X FUNCTIONAL TOTAL ENERGY']}
pv0['PBE0 TOTAL ENERGY'] = {'func': sum, 'args': ['PBE0 FUNCTIONAL TOTAL ENERGY']}
pv0['PBE TOTAL ENERGY'] = {'func': sum, 'args': ['PBE FUNCTIONAL TOTAL ENERGY']}
pv0['HF TOTAL ENERGY'] = {'func': sum, 'args': ['SCF TOTAL ENERGY']}
pv0['CCSD(T) CORRELATION ENERGY'] = {'func': difference, 'args': ['CCSD(T) FNO CORRELATION ENERGY', 'FNO CORRECTION ENERGY']}
pv0['CCSD(T) TOTAL ENERGY'] = {'func': sum, 'args': ['HF TOTAL ENERGY', 'CCSD(T) CORRELATION ENERGY']}
pv0['CCSDT(Q) TOTAL ENERGY'] = {'func': sum, 'args': ['HF TOTAL ENERGY', 'CCSDT(Q) CORRELATION ENERGY']}
pv0['VV10 TOTAL ENERGY'] = {'func': sum, 'args': ['VV10 FUNCTIONAL TOTAL ENERGY']}
pv0['LC-VV10 TOTAL ENERGY'] = {'func': sum, 'args': ['LC-VV10 FUNCTIONAL TOTAL ENERGY']}
pv0['M08-HX TOTAL ENERGY'] = {'func': sum, 'args': ['M08-HX FUNCTIONAL TOTAL ENERGY']}
pv0['M08-SO TOTAL ENERGY'] = {'func': sum, 'args': ['M08-SO FUNCTIONAL TOTAL ENERGY']}
pv0['M11 TOTAL ENERGY'] = {'func': sum, 'args': ['M11 FUNCTIONAL TOTAL ENERGY']}
pv0['M11L TOTAL ENERGY'] = {'func': sum, 'args': ['M11L FUNCTIONAL TOTAL ENERGY']}
pv0['MP2 TOTAL ENERGY'] = {'func': sum, 'args':['HF TOTAL ENERGY', 'MP2 CORRELATION ENERGY']}
pv0['CCSD TOTAL ENERGY'] = {'func': sum, 'args':['HF TOTAL ENERGY', 'CCSD CORRELATION ENERGY']}

for pvar, action in pv0.iteritems():
    try:
        df = append_result_of_func_with_funcargs_to_master_DataFrame_atlevel_with_label(
            master=df, label=pvar, atlevel=lvl, func=action['func'], funcargs=action['args'])
        print """Built df pv0: %s""" % (pvar)
        #print df.xs(pvar, level=lvl).head(2)
    except KeyError, e:
        print """Error building df pv0: empty index '%s' because missing %s""" % (pvar, e)

# <<< miscellaneous pre-computing >>>
try:
    df.xs('HF-CABS TOTAL ENERGY', level='psivar')
    df.xs('MP2-F12 TOTAL ENERGY', level='psivar')
except KeyError, e:
    print 'Not handled: DF_OMEGA', e
else:
    ratio_HFCABS_MP2F12 = ie(df.xs('HF-CABS TOTAL ENERGY', level='psivar')) / ie(df.xs('MP2-F12 TOTAL ENERGY', level='psivar'))
    #ratio_HFCABS_MP2F12_unCP = ie_uncp(df.xs('HF-CABS TOTAL ENERGY', level='psivar')) / ie(df.xs('MP2-F12 TOTAL ENERGY', level='psivar'))
    print ratio_HFCABS_MP2F12
    df_omega = pd.concat({
        'adz': omega([-1, 4, ratio_HFCABS_MP2F12['adz']]),
        'atz': omega([0.4, 0.6, ratio_HFCABS_MP2F12['atz']]),
        })
    df_omega = pd.DataFrame(df_omega, columns=['dimer'])
    df_omega['psivar'] = 'DW-CCSD(T)-F12 OMEGA'
    df_omega['monoA'] = df_omega['dimer']
    df_omega['monoB'] = df_omega['dimer']
    df_omega.set_index('psivar', append=True, inplace=True)
    df_omega_uncp = df_omega.copy()  # utterly wrong! TODO
    df_omega2 = pd.concat([df_omega, df_omega_uncp], keys=['CP','unCP'], axis=1)
    df_omega = df_omega2.copy()
    df_omega = df_omega.reorder_levels([0, 3, 1, 2])
    df_omega.index.names = ['bstrt', 'psivar', 'meta', 'rxn']
    df_omega['unCP'] = np.nan
    df = df.append(df_omega, verify_integrity=True)

try:
    df.xs('nobas', level='bstrt')
except KeyError, e:
    print 'Not handled: DASH-D', e
else:
    df_nobas = pd.concat({tup[0]: df.xs('nobas', level='bstrt') for tup in df.index.values if tup[0] != 'nobas'})
    df_nobas.index.names = ['bstrt', 'psivar', 'rxn']
    df = df.append(df_nobas, verify_integrity=True)

# <<< append to main DataFrame computable method quantities >>>

lvl = 'psivar'
pv1 = collections.OrderedDict()
pv1['DW-CCSD(T**)-F12 CORRELATION ENERGY'] = {'func': dispersion_weighting, 'args': ['DW-CCSD(T)-F12 OMEGA', 'CCSD(T**)-F12A CORRELATION ENERGY', 'CCSD(T**)-F12B CORRELATION ENERGY']}
pv1['DW-CCSD(T**)-F12 TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'DW-CCSD(T**)-F12 CORRELATION ENERGY']}
pv1['B3LYP-D2 TOTAL ENERGY'] = {'func': sum, 'args': ['B3LYP FUNCTIONAL TOTAL ENERGY', 'B3LYP-D2 DISPERSION CORRECTION ENERGY']}
pv1['B3LYP-D3 TOTAL ENERGY'] = {'func': sum, 'args': ['B3LYP FUNCTIONAL TOTAL ENERGY', 'B3LYP-D3 DISPERSION CORRECTION ENERGY']}
pv1['B3LYP-D3BJ TOTAL ENERGY'] = {'func': sum, 'args': ['B3LYP FUNCTIONAL TOTAL ENERGY', 'B3LYP-D3BJ DISPERSION CORRECTION ENERGY']}
pv1['B3LYP-XDM TOTAL ENERGY'] = {'func': sum, 'args': ['B3LYP FUNCTIONAL TOTAL ENERGY', 'B3LYP-XDM DISPERSION CORRECTION ENERGY']}
pv1['B2PLYP-D2 TOTAL ENERGY'] = {'func': sum, 'args': ['B2PLYP TOTAL ENERGY', 'B2PLYP-D2 DISPERSION CORRECTION ENERGY']}
pv1['B2PLYP-D3 TOTAL ENERGY'] = {'func': sum, 'args': ['B2PLYP TOTAL ENERGY', 'B2PLYP-D3 DISPERSION CORRECTION ENERGY']}
pv1['B2PLYP-D3BJ TOTAL ENERGY'] = {'func': sum, 'args': ['B2PLYP TOTAL ENERGY', 'B2PLYP-D3BJ DISPERSION CORRECTION ENERGY']}
pv1['B970-D2 TOTAL ENERGY'] = {'func': sum, 'args': ['B970 FUNCTIONAL TOTAL ENERGY', 'B970-D2 DISPERSION CORRECTION ENERGY']}
pv1['B97-D2 TOTAL ENERGY'] = {'func': sum, 'args': ['B97 FUNCTIONAL TOTAL ENERGY', 'B97-D2 DISPERSION CORRECTION ENERGY']}
pv1['B97-D3 TOTAL ENERGY'] = {'func': sum, 'args': ['B97 FUNCTIONAL TOTAL ENERGY', 'B97-D3 DISPERSION CORRECTION ENERGY']}
pv1['B97-D3BJ TOTAL ENERGY'] = {'func': sum, 'args': ['B97 FUNCTIONAL TOTAL ENERGY', 'B97-D3BJ DISPERSION CORRECTION ENERGY']}
pv1['BP86-D2 TOTAL ENERGY'] = {'func': sum, 'args': ['BP86 FUNCTIONAL TOTAL ENERGY', 'BP86-D2 DISPERSION CORRECTION ENERGY']}
pv1['BP86-D3 TOTAL ENERGY'] = {'func': sum, 'args': ['BP86 FUNCTIONAL TOTAL ENERGY', 'BP86-D3 DISPERSION CORRECTION ENERGY']}
pv1['BP86-D3BJ TOTAL ENERGY'] = {'func': sum, 'args': ['BP86 FUNCTIONAL TOTAL ENERGY', 'BP86-D3BJ DISPERSION CORRECTION ENERGY']}
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
for pvar, action in pv1.iteritems():
    try:
        df = append_result_of_func_with_funcargs_to_master_DataFrame_atlevel_with_label(
            master=df, label=pvar, atlevel=lvl, func=action['func'], funcargs=action['args'])
        #print df.xs(pvar, level=lvl).head(2)
    except KeyError, e:
        print """Error building df pv0: empty index '%s' because missing %s""" % (pvar, e)

#print {(tup[0], tup[1]) for tup in df.index.values if 'TOTAL' in tup[1]}
#print df.xs('adz', level='bstrt').xs('B3LYP TOTAL ENERGY', level='psivar')
#print df.xs('adz', level='bstrt').xs('B3LYP-D2 DISPERSION CORRECTION ENERGY', level='psivar')

# <<< define extra Method and BasisSet objects >>>

# Note dict key must match object name; i.e., first two strings on line must be the same
#   Note also that methods are all caps, bases all lowercase

#bases['hill1_adtz'] = BasisSet('hill1_adtz', build=[['hillcc_adtz'], ['atz', 'hillcc_adtz']])  # TODO should have None or non-xtpl first element?
#bases['hill2_dtzf12'] = BasisSet('hill2_dtzf12', build=[None, None, ['tzf12', 'hillcc_dtzf12', 'hillt_dtzf12']])
#methods['CCSDTNSAF12'] = Method('CCSDTNSAF12', fullname='CCSD(T)-F12a')

# <<< append to main DataFrame computable basis treatment quantities >>>

lvl = 'bstrt'
pv2 = collections.OrderedDict()
pv2['adtz'] = {'func': xtpl_power, 'args': [3.0, 3, 'atz', 'adz']}
pv2['atqz'] = {'func': xtpl_power, 'args': [3.0, 4, 'aqz', 'atz']}
pv2['aq5z'] = {'func': xtpl_power, 'args': [3.0, 5, 'a5z', 'aqz']}
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
        df = append_result_of_func_with_funcargs_to_master_DataFrame_atlevel_with_label(
            master=df, label=pvar, atlevel=lvl, func=action['func'], funcargs=action['args'])
        print """Built df pv2: %s""" % (pvar)
        #print df.xs(pvar, level=lvl).head(2)
    except KeyError, e:
        print """Error building df pv2: empty index '%s' because missing %s""" % (pvar, e)

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
    elif 'MP2' in stub:
        return 2
    else:
        return 3

def generic_mtd(Nstage, stub):
    if Nstage == 1:
        return ['%s TOTAL ENERGY' % (stub)]
    elif Nstage == 2:
        return ['%s TOTAL ENERGY' % ('HF-CABS' if 'F12' in stub else 'SCF'), '%s CORRELATION ENERGY' % (stub)]
    elif Nstage == 3:
        #return ['%s TOTAL ENERGY' % ('HF-CABS' if 'F12' in stub else 'SCF'), '%s CORRELATION ENERGY' % (stub), '%s CC CORRECTION ENERGY' % (stub)]
        return ['%s TOTAL ENERGY' % ('HF-CABS' if 'F12' in stub else 'SCF'), '%s CORRELATION ENERGY' % ('MP2-F12' if 'F12' in stub else 'MP2'), 
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
    cpdict = {'CP': ie, 'unCP': ie_uncp, 'ave': ie_ave}
    func = cpdict[cpmode]
    #func = ie if cpmode == 'CP' else ie_uncp
    if baslist is None:
        raise KeyError  # TODO a more specific message that mtd/bas don't mix wouldn't hurt
        #print '\n <<<', methods[method].fullname, '/', bases[basis].fullname, '>>>'
        #print 'stages:', 'M:', compute_max_mtd(method), 'B:', compute_max_bas(basis), 'U:', Nstage
    #if method == 'B2PLYPD3' and option == 'nfc' and cpmode == 'CP' and basis == 'adz':
    #if basis == 'atqzadz':
        #print 'pcss:', mtdlist
        #print 'bass:', baslist
        #print 'opts:', optlist
        #for pcs, bas, opt in zip(mtdlist, baslist, optlist):
            #print df.loc[bas].loc[pcs].loc[opt].loc['A24-4']
    return func(sum([df.loc[bas].loc[pcs].loc[opt] for pcs, bas, opt in zip(mtdlist, baslist, optlist)]))


# <<< assemble all model chemistries into columns of new DataFrame >>>

rxns = ['%s-%s' % (dbse, rxn) for rxn in dbobj.hrxn.keys()]
mine = pd.DataFrame({}, index=rxns)
mine.index.names = ['rxn']

#print df.xs('aqz', level='bstrt').xs('HF TOTAL ENERGY', level='psivar').xs('A24-4', level='rxn')
#print df.xs('atqz', level='bstrt').xs('MP2 CORRELATION ENERGY', level='psivar').xs('A24-4', level='rxn')
#print df.xs('adz', level='bstrt').xs('CCSD(T) CC CORRECTION ENERGY', level='psivar').xs('A24-4', level='rxn')

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
    mtds = ['HF', 'CCSDT', 'CCSDTQ']
    bass = ['adz', 'atz', 'hadz', 'jadz', 'atz', 'adtz', 'aq5z']
    opts = ['', 'full', 'fno1e5', 'fno1e4', 'fno5e5', 'fno1e6', 'fno1e4-mrcc']
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

for cpm in cpmd:
    for mtd in mtds:
        for bas in bass:
            for opt in opts:
                mc = '-'.join([mtd, cpm, bas]) if opt == '' else '-'.join([mtd, opt, cpm, bas])
                try:
                    mine[mc] = build(mtd, opt, cpm, bas)
                    print """Built model chemistry %s""" % (mc)
                except KeyError, e:
                    print """Error building model chemistry: empty index '%s' because missing %s""" % (mc, e)
                    pass

def threadtheframe(modelchem, xlimit=4.0):
    dbdat = []
    for rxn in dbobj.hrxn.keys():
        data = []
        for mc in modelchem:
            value = h2kc * mine[mc]['%s-%s' % (dbse, rxn)]
            #if pd.isnull(value):
            #    data.append(None)
            #else:
            #    data.append(value)
            data.append(None if pd.isnull(value) else value)
        dbdat.append({'sys': str(rxn), 'color': dbobj.hrxn[rxn].color, 'data': data})
    mpl.thread(dbdat, modelchem, color='sapt', xlimit=xlimit)

if project == 'f12dilabio':
    mine['CCSDTNSBF12-CP-hill2_adtz'] = build_from_lists(['HF-CABS TOTAL ENERGY', 'CCSD-F12B CORRELATION ENERGY', '(T)-F12AB CORRECTION ENERGY'], ['atz', 'hillcc_adtz', 'hillt_adtz'])
    
if project == 'parenq':
    # TODO these aren't actually interaction energies!
    mine['CCSDTQ-corl-CP-adz'] = build_from_lists(['CCSDT(Q) CORRELATION ENERGY'], ['adz'], ['fno1e4'])
    mine['CCSDT-corl-CP-adz'] = build_from_lists(['CCSD(T) CORRELATION ENERGY'], ['adz'], ['fno1e4'])
    mine['CCSD-corl-CP-adz'] = build_from_lists(['CCSD CORRELATION ENERGY'], ['adz'], ['fno1e4-mrcc'])
    mine['MP2-corl-CP-adz'] = build_from_lists(['MP2 CORRELATION ENERGY'], ['adz'])
    mine['CCSDTQ-corl-CP-atz'] = build_from_lists(['CCSDT(Q) CORRELATION ENERGY'], ['atz'], ['fno1e4'])
    mine['CCSDT-corl-CP-atz'] = build_from_lists(['CCSD(T) CORRELATION ENERGY'], ['atz'], ['fno1e4'])
    mine['CCSD-corl-CP-atz'] = build_from_lists(['CCSD CORRELATION ENERGY'], ['atz'], ['fno1e4'])
    mine['MP2-corl-CP-atz'] = build_from_lists(['MP2 CORRELATION ENERGY'], ['atz'])
    
#print mine.columns
#print mine['CCSDTQ-corl-CP-adz']
#print mine['CCSDT-corl-CP-adz']
#print mine['CCSD-corl-CP-adz']
#print mine['MP2-corl-CP-adz']
    mine['convdel-adz'] = mine['CCSDT-corl-CP-adz'] - mine['MP2-corl-CP-adz']
    mine['altdel-adz'] = mine['CCSDTQ-corl-CP-adz'] - mine['CCSD-corl-CP-adz']
    minelist = ['convdel-adz', 'altdel-adz', 'CCSDTQ-corl-CP-adz', 'CCSDT-corl-CP-adz', 'CCSD-corl-CP-adz', 'MP2-corl-CP-adz']
    mine['convdel-atz'] = mine['CCSDT-corl-CP-atz'] - mine['MP2-corl-CP-atz']
    mine['altdel-atz'] = mine['CCSDTQ-corl-CP-atz'] - mine['CCSD-corl-CP-atz']
    minelistatz = ['convdel-atz', 'altdel-atz', 'CCSDTQ-corl-CP-atz', 'CCSDT-corl-CP-atz', 'CCSD-corl-CP-atz', 'MP2-corl-CP-atz']

    #import mpl
    #import matplotlib.pyplot as plt
    #threadtheframe(minelist)
    #threadtheframe(minelistatz)
    #threadtheframe(['CCSDTQ-fno1e4-CP-adz', 'CCSDTQ-fno1e4-CP-atz', 'CCSDTQ-fno1e4-CP-adtz', 'CCSDTQ-fno1e5-CP-adtz', 'CCSDTQ-fno1e5-CP-atz', 'CCSDTQ-fno1e5-CP-adz', 'CCSDT-CP-aq5z'], xlimit=7.0)




# <<< test cases >>>

from qcdb.pandas_test import test_pandas
test_pandas(h2kc, project, mine)

# <<< write qcdb data loader >>>

f1 = open('%s_%s.py' % (dbse, project), 'w')
f1.write('\ndef load_%s(dbinstance):\n\n' % (project))

for mc in mine.columns:
    lmc = mc.split('-')  # TODO could be done better
    method = lmc[0]
    bsse = '_'.join(lmc[1:-1])
    basis = lmc[-1]
    for rxn in dbobj.hrxn.keys():
        value = h2kc * mine[mc]['%s-%s' % (dbse, rxn)]
        if pd.isnull(value):
            pass
        else:
            f1.write("""    dbinstance.add_ReactionDatum(dbse='%s', rxn=%s, method='%s', mode='%s', basis='%s', value=%.4f)\n""" %
                (dbse, repr(rxn), method, bsse, basis, value))

f1.close()

            
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
